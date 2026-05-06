/*
 * Dispatch layer between htslib's sam_/hts_ API and libisam.
 *
 * Strategy: tagged pointers. Real htslib pointers (htsFile/hts_idx_t/hts_itr_t)
 * come from malloc and are >=8-byte aligned, so bit 0 is always clear. We OR
 * in 1 to mark our own pointers and untag at the dispatch boundary. No side
 * table, no mutex, no cost on the htslib path.
 *
 * The file type is decided once at hts_open / sam_open time by looking at
 * the path suffix; the resulting tag carries through to hdr / index / iterator
 * pointers handed back to the caller.
 */
#define HTSISAM_NO_DROPIN
#include "../include/isam/isam.h"
#include "../include/htsisam/sam.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>

#define ISAM_TAG     ((uintptr_t)1)
#define IS_ISAM(p)   ((p) && ((uintptr_t)(p) & ISAM_TAG))
#define TAG(p)       ((void *)((uintptr_t)(p) | ISAM_TAG))
#define UNTAG(p)     ((void *)((uintptr_t)(p) & ~ISAM_TAG))

/* The wrappers cast tagged isam objects through generic htslib types so the
 * caller can store them in samFile / hts_idx_t / hts_itr_t pointers without
 * knowing the difference. */

int htsisam_is_isam_path(const char *path) {
    if (!path) return 0;
    size_t n = strlen(path);
    if (n < 5) return 0;
    return strcasecmp(path + n - 5, ".isam") == 0;
}

/* ---- open / close ----------------------------------------------------- */

static htsFile *open_dispatch(const char *fn, const char *mode) {
    if (htsisam_is_isam_path(fn)) {
        if (!mode || mode[0] != 'r') {
            /* writing into a .isam isn't a sam_write1 operation; it's an
             * import. Refuse rather than silently produce wrong output. */
            fprintf(stderr,
                    "htsisam: write mode '%s' is not supported for .isam "
                    "(use `isamtools import`)\n", mode);
            return NULL;
        }
        isamFile *f = isam_open(fn, "r");
        if (!f) return NULL;
        return (htsFile *)TAG(f);
    }
    return hts_open(fn, mode);
}

htsFile *htsisam_hts_open(const char *fn, const char *mode) {
    return open_dispatch(fn, mode);
}

htsFile *htsisam_sam_open(const char *fn, const char *mode) {
    return open_dispatch(fn, mode);
}

htsFile *htsisam_hts_open_format(const char *fn, const char *mode,
                                 const htsFormat *fmt) {
    if (htsisam_is_isam_path(fn))
        return open_dispatch(fn, mode);
    return hts_open_format(fn, mode, fmt);
}

int htsisam_sam_close(htsFile *fp) {
    if (IS_ISAM(fp)) return isam_close((isamFile *)UNTAG(fp));
    return sam_close(fp);
}

int htsisam_hts_close(htsFile *fp) {
    if (IS_ISAM(fp)) return isam_close((isamFile *)UNTAG(fp));
    return hts_close(fp);
}

/* ---- header ----------------------------------------------------------- */

sam_hdr_t *htsisam_sam_hdr_read(htsFile *fp) {
    if (IS_ISAM(fp)) return isam_hdr_read((isamFile *)UNTAG(fp));
    return sam_hdr_read(fp);
}

int htsisam_sam_hdr_write(htsFile *fp, const sam_hdr_t *hdr) {
    if (IS_ISAM(fp)) {
        fprintf(stderr, "htsisam: sam_hdr_write to .isam is not supported\n");
        return -1;
    }
    return sam_hdr_write(fp, hdr);
}

/* ---- record IO -------------------------------------------------------- */

int htsisam_sam_read1(htsFile *fp, sam_hdr_t *hdr, bam1_t *b) {
    if (IS_ISAM(fp)) return isam_read1((isamFile *)UNTAG(fp), hdr, b);
    return sam_read1(fp, hdr, b);
}

int htsisam_sam_write1(htsFile *fp, const sam_hdr_t *hdr, const bam1_t *b) {
    if (IS_ISAM(fp)) {
        fprintf(stderr, "htsisam: sam_write1 to .isam is not supported\n");
        return -1;
    }
    return sam_write1(fp, hdr, b);
}

/* ---- index -----------------------------------------------------------
 *
 * For .isam files the "index" is the file itself; we open a separate read
 * handle on the same path. For BAM/CRAM we delegate to htslib so the
 * .bai/.csi/.crai sidecars are loaded normally.
 */

hts_idx_t *htsisam_sam_index_load(htsFile *fp, const char *fn) {
    if (IS_ISAM(fp)) {
        isamIndex *i = isam_index_load(fn);
        return i ? (hts_idx_t *)TAG(i) : NULL;
    }
    return sam_index_load(fp, fn);
}

hts_idx_t *htsisam_sam_index_load2(htsFile *fp, const char *fn,
                                   const char *fnidx) {
    if (IS_ISAM(fp)) {
        isamIndex *i = isam_index_load(fnidx ? fnidx : fn);
        return i ? (hts_idx_t *)TAG(i) : NULL;
    }
    return sam_index_load2(fp, fn, fnidx);
}

hts_idx_t *htsisam_sam_index_load3(htsFile *fp, const char *fn,
                                   const char *fnidx, int flags) {
    if (IS_ISAM(fp)) {
        isamIndex *i = isam_index_load(fnidx ? fnidx : fn);
        return i ? (hts_idx_t *)TAG(i) : NULL;
    }
    return sam_index_load3(fp, fn, fnidx, flags);
}

void htsisam_hts_idx_destroy(hts_idx_t *idx) {
    if (IS_ISAM(idx)) { isam_index_destroy((isamIndex *)UNTAG(idx)); return; }
    hts_idx_destroy(idx);
}

/* ---- iterators -------------------------------------------------------- */

hts_itr_t *htsisam_sam_itr_querys(const hts_idx_t *idx, sam_hdr_t *hdr,
                                  const char *region) {
    if (IS_ISAM(idx)) {
        isamItr *it = isam_itr_querys((isamIndex *)UNTAG(idx), hdr, region);
        return it ? (hts_itr_t *)TAG(it) : NULL;
    }
    return sam_itr_querys(idx, hdr, region);
}

hts_itr_t *htsisam_sam_itr_queryi(const hts_idx_t *idx, int tid,
                                  hts_pos_t beg, hts_pos_t end) {
    if (IS_ISAM(idx)) {
        isamItr *it = isam_itr_queryi((isamIndex *)UNTAG(idx), tid, beg, end);
        return it ? (hts_itr_t *)TAG(it) : NULL;
    }
    return sam_itr_queryi(idx, tid, beg, end);
}

int htsisam_sam_itr_next(htsFile *fp, hts_itr_t *itr, bam1_t *b) {
    if (IS_ISAM(itr)) {
        if (!IS_ISAM(fp)) {
            fprintf(stderr,
                    "htsisam: itr is from a .isam but fp is plain BAM/CRAM\n");
            return -2;
        }
        return isam_itr_next((isamFile *)UNTAG(fp),
                             (isamItr *)UNTAG(itr), b);
    }
    /* Replicate htslib's static-inline body. We can't call the original
     * sam_itr_next because the header redirected it to ourselves. */
    if (!fp->is_bgzf && !fp->is_cram) return -2;
    if (!itr) return -2;
    if (itr->multi)
        return hts_itr_multi_next(fp, itr, b);
    return hts_itr_next(fp->is_bgzf ? fp->fp.bgzf : NULL, itr, b, fp);
}

void htsisam_hts_itr_destroy(hts_itr_t *itr) {
    if (IS_ISAM(itr)) { isam_itr_destroy((isamItr *)UNTAG(itr)); return; }
    hts_itr_destroy(itr);
}
