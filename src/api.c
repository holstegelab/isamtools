#include "../include/isam/isam.h"
#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"
#include "aux_tags.h"

#include <htslib/hts.h>

struct isamFile {
    view_t       *v;
    lmdb_view_iter_t *seq_lmdb;
    bam1_t       *canon;
    payload_cache_t *cache;
    isam_tag_set_t shared_tags;
    isam_tag_set_t dropped_tags;
    int           errcode;
};

struct isamIndex {
    char    *path;
    lmdb_backend_t *lmdb;
};

struct isamItr {
    lmdb_view_iter_t *lmdb_it;
    int           done;
    int           filter_region;
    int           tid;
    int64_t       beg;
    int64_t       end;
};

static int record_to_bam(isamFile *fp,
                         const lmdb_view_record_t *rec,
                         bam1_t *b) {
    store_record_t srec;
    if (!store_lookup_read_id(fp->v->store, rec->store_read_id, &srec)) {
        fp->errcode = ISAM_ERR_MISSING_QNAME;
        return ISAM_ERR_MISSING_QNAME;
    }
    int frc = payload_cache_fetch(fp->cache, fp->v->store, &srec, fp->canon);
    if (frc != 0) {
        store_record_clear(&srec);
        fp->errcode = frc;
        return frc;
    }

    isam_rebuild_bam(b, fp->canon, rec->flag, rec->tid,
                     rec->pos0, rec->mapq, rec->n_cigar, rec->cigar,
                     rec->mate_tid, rec->mate_pos0, rec->tlen,
                     rec->aux, rec->aux_len);
    lmdb_read_aux_apply(fp->v->lmdb, rec->store_read_id, b, 0);
    lmdb_read_tags_apply(fp->v->lmdb, rec->store_read_id, b, 0);
    isam_aux_tags_append_missing_from(fp->canon, b,
                                      &fp->shared_tags,
                                      &fp->dropped_tags);
    lmdb_read_tags_apply(fp->v->store->lmdb, rec->store_read_id, b, 0);
    isam_aux_tags_remove(b, &fp->dropped_tags);
    store_record_clear(&srec);
    /* Match htslib: positive on success (b->l_data is the record's serialized
     * size, always > 0 for any non-empty record). Consumers that test `> 0`
     * (a common idiom) work without modification. */
    return b->l_data > 0 ? b->l_data : 1;
}

isamFile *isam_open(const char *path, const char *mode) {
    return isam_open_ext(path, mode, NULL);
}

isamFile *isam_open_ext(const char *path, const char *mode,
                        const isamOptions *opts) {
    if (!mode || mode[0] != 'r') return NULL;
    isamFile *fp = xcalloc(1, sizeof *fp);
    fp->v = view_open_read(path);
    view_link_store(fp->v,
                    opts ? opts->store_path_override       : NULL,
                    opts ? opts->store_index_path_override : NULL);
    if (opts && opts->reference_fasta &&
        hts_set_fai_filename(fp->v->store->sam, opts->reference_fasta) != 0) {
        view_close(fp->v);
        free(fp);
        return NULL;
    }
    if (opts && opts->n_threads > 1)
        hts_set_threads(fp->v->store->sam, opts->n_threads);
    fp->canon = bam_init1();
    fp->cache = payload_cache_create(opts ? opts->payload_cache_bytes : 0);
    isam_tag_set_init(&fp->shared_tags);
    isam_tag_set_init(&fp->dropped_tags);
    char *shared_csv = view_meta_get(fp->v, "shared_tags_snapshot");
    char *dropped_csv = view_meta_get(fp->v, "dropped_tags");
    if (isam_tag_set_parse_csv(&fp->shared_tags, shared_csv) != 0 ||
        isam_tag_set_parse_csv(&fp->dropped_tags, dropped_csv) != 0) {
        free(shared_csv);
        free(dropped_csv);
        isam_close(fp);
        return NULL;
    }
    free(shared_csv);
    free(dropped_csv);
    return fp;
}

int isam_close(isamFile *fp) {
    if (!fp) return 0;
    if (fp->seq_lmdb) lmdb_view_iter_close(fp->seq_lmdb);
    if (fp->canon)    bam_destroy1(fp->canon);
    if (fp->cache)    payload_cache_destroy(fp->cache);
    isam_tag_set_destroy(&fp->shared_tags);
    isam_tag_set_destroy(&fp->dropped_tags);
    if (fp->v)        view_close(fp->v);
    free(fp);
    return 0;
}

sam_hdr_t *isam_hdr_read(isamFile *fp) {
    return sam_hdr_dup(fp->v->hdr);
}

int isam_hdr_write(isamFile *fp, const sam_hdr_t *hdr) {
    (void)fp; (void)hdr;
    return ISAM_ERR_UNSUPPORTED;
}

int isam_read1(isamFile *fp, sam_hdr_t *hdr, bam1_t *b) {
    (void)hdr;
    if (!fp->seq_lmdb) {
        int64_t has_coord = 0;
        view_meta_get_i64(fp->v, "indexed_coord", &has_coord);
        fp->seq_lmdb = lmdb_view_iter_open(fp->v->lmdb,
            has_coord ? LMDB_VIEW_ITER_COORD : LMDB_VIEW_ITER_IMPORT,
            -1, 0, 0);
    }
    lmdb_view_record_t rec;
    if (!lmdb_view_iter_next(fp->seq_lmdb, &rec)) return ISAM_EOF;
    int rc = record_to_bam(fp, &rec, b);
    lmdb_view_record_clear(&rec);
    return rc;
}

isamIndex *isam_index_load(const char *path) {
    isamIndex *idx = xcalloc(1, sizeof *idx);
    idx->path = xstrdup(path);
    idx->lmdb = lmdb_backend_open(path, /*readonly*/1, 0);
    return idx;
}

void isam_index_destroy(isamIndex *idx) {
    if (!idx) return;
    if (idx->lmdb) lmdb_backend_close(idx->lmdb);
    free(idx->path);
    free(idx);
}

isamItr *isam_itr_queryi(isamIndex *idx, int tid, hts_pos_t beg, hts_pos_t end) {
    int64_t has_iv = 0;
    lmdb_backend_meta_get_i64(idx->lmdb, "indexed_interval", &has_iv);
    if (!has_iv) return NULL;

    isamItr *it = xcalloc(1, sizeof *it);
    it->lmdb_it = lmdb_view_iter_open(idx->lmdb, LMDB_VIEW_ITER_REGION,
                                      tid, (int64_t)beg, (int64_t)end);
    return it;
}

isamItr *isam_itr_querys(isamIndex *idx, sam_hdr_t *hdr, const char *region) {
    int tid = -1;
    hts_pos_t beg = 0, end = HTS_POS_MAX;
    const char *tail = sam_parse_region(hdr, region, &tid, &beg, &end, 0);
    if (!tail || tid < 0) return NULL;
    return isam_itr_queryi(idx, tid, beg, end);
}

isamItr *isam_itr_query_name(isamIndex *idx, const char *qname) {
    isamItr *it = xcalloc(1, sizeof *it);
    int64_t store_read_id = 0;
    if (lmdb_view_lookup_qname(idx->lmdb, qname_hash(qname), qname,
                               &store_read_id)) {
        it->lmdb_it = lmdb_view_iter_open_read_id(idx->lmdb, store_read_id);
    } else {
        it->done = 1;
    }
    return it;
}

isamItr *isam_itr_query_subset(isamIndex *idx, const char *subset_name) {
    if (!subset_name) return NULL;
    int64_t subset_id = 0;
    if (!lmdb_subset_lookup(idx->lmdb, subset_name, &subset_id)) return NULL;
    isamItr *it = xcalloc(1, sizeof *it);
    it->lmdb_it = lmdb_view_iter_open_subset(idx->lmdb, subset_id);
    return it;
}

isamItr *isam_itr_query_region_subset(isamIndex *idx,
                                      sam_hdr_t *hdr,
                                      const char *region,
                                      const char *subset_name) {
    if (!subset_name) return isam_itr_querys(idx, hdr, region);
    int tid = -1;
    hts_pos_t beg = 0, end = HTS_POS_MAX;
    const char *tail = sam_parse_region(hdr, region, &tid, &beg, &end, 0);
    if (!tail || tid < 0) return NULL;

    isamItr *it = isam_itr_query_subset(idx, subset_name);
    if (!it) return NULL;
    it->filter_region = 1;
    it->tid = tid;
    it->beg = (int64_t)beg;
    it->end = (int64_t)end;
    return it;
}

int isam_itr_next(isamFile *fp, isamItr *itr, bam1_t *b) {
    if (!itr || itr->done) return ISAM_EOF;
    lmdb_view_record_t rec;
    while (lmdb_view_iter_next(itr->lmdb_it, &rec)) {
        if (itr->filter_region &&
            !(rec.tid == itr->tid && rec.pos0 < itr->end && rec.end0 > itr->beg)) {
            lmdb_view_record_clear(&rec);
            continue;
        }
        int rc = record_to_bam(fp, &rec, b);
        lmdb_view_record_clear(&rec);
        return rc;
    }
    itr->done = 1;
    return ISAM_EOF;
}

void isam_itr_destroy(isamItr *itr) {
    if (!itr) return;
    if (itr->lmdb_it) lmdb_view_iter_close(itr->lmdb_it);
    free(itr);
}

const char *isam_strerror(int errcode) {
    switch (errcode) {
    case ISAM_OK:                  return "ok";
    case ISAM_EOF:                 return "end of stream";
    case ISAM_ERR_IO:              return "i/o error";
    case ISAM_ERR_BACKEND:         return "backend error";
    case ISAM_ERR_FORMAT:          return "format error";
    case ISAM_ERR_MISSING_STORE:   return "canonical store not available";
    case ISAM_ERR_MISSING_QNAME:   return "qname missing from store";
    case ISAM_ERR_STALE_INDEX:     return "store/view index is stale";
    case ISAM_ERR_NO_INDEX:        return "required index is missing";
    case ISAM_ERR_UNSUPPORTED:     return "operation not supported";
    default:                       return "unknown error";
    }
}

int isam_errno(isamFile *fp) {
    return fp ? fp->errcode : 0;
}
