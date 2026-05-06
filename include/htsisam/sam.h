/*
 * htsisam: drop-in compatibility layer letting htslib-using code transparently
 * read .isam files alongside regular SAM/BAM/CRAM. Include after (or in place
 * of) <htslib/sam.h>; the macros below redirect a small set of htslib entry
 * points to dispatching wrappers that detect the file type at runtime.
 *
 * Suggested integration: add `-include htsisam/sam.h` to your CFLAGS so every
 * TU sees the redirects without source edits, then swap `-lhts` for
 * `-lisam -lhts` at link time.
 */
#ifndef HTSISAM_SAM_H
#define HTSISAM_SAM_H

#include <htslib/hts.h>
#include <htslib/sam.h>

#if defined(_WIN32)
#  define HTSISAM_EXPORT __declspec(dllexport)
#else
#  define HTSISAM_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

HTSISAM_EXPORT htsFile  *htsisam_hts_open       (const char *fn, const char *mode);
HTSISAM_EXPORT htsFile  *htsisam_hts_open_format(const char *fn, const char *mode,
                                                 const htsFormat *fmt);
HTSISAM_EXPORT htsFile  *htsisam_sam_open       (const char *fn, const char *mode);
HTSISAM_EXPORT int       htsisam_sam_close      (htsFile *fp);
HTSISAM_EXPORT int       htsisam_hts_close      (htsFile *fp);

HTSISAM_EXPORT sam_hdr_t *htsisam_sam_hdr_read  (htsFile *fp);
HTSISAM_EXPORT int        htsisam_sam_hdr_write (htsFile *fp, const sam_hdr_t *hdr);

HTSISAM_EXPORT int        htsisam_sam_read1     (htsFile *fp, sam_hdr_t *hdr, bam1_t *b);
HTSISAM_EXPORT int        htsisam_sam_write1    (htsFile *fp, const sam_hdr_t *hdr,
                                                 const bam1_t *b);

HTSISAM_EXPORT hts_idx_t *htsisam_sam_index_load (htsFile *fp, const char *fn);
HTSISAM_EXPORT hts_idx_t *htsisam_sam_index_load2(htsFile *fp, const char *fn,
                                                  const char *fnidx);
HTSISAM_EXPORT hts_idx_t *htsisam_sam_index_load3(htsFile *fp, const char *fn,
                                                  const char *fnidx, int flags);
HTSISAM_EXPORT void       htsisam_hts_idx_destroy(hts_idx_t *idx);

HTSISAM_EXPORT hts_itr_t *htsisam_sam_itr_querys (const hts_idx_t *idx, sam_hdr_t *hdr,
                                                  const char *region);
HTSISAM_EXPORT hts_itr_t *htsisam_sam_itr_queryi (const hts_idx_t *idx, int tid,
                                                  hts_pos_t beg, hts_pos_t end);
HTSISAM_EXPORT int        htsisam_sam_itr_next   (htsFile *fp, hts_itr_t *itr, bam1_t *b);
HTSISAM_EXPORT void       htsisam_hts_itr_destroy(hts_itr_t *itr);

/* Returns 1 if the path resolves to a .isam file (suffix or magic). */
HTSISAM_EXPORT int        htsisam_is_isam_path(const char *path);

#ifdef __cplusplus
}
#endif

/* Drop-in macros. Define HTSISAM_NO_DROPIN before including this header to
 * keep the htslib symbols pristine and call htsisam_* explicitly. */
#ifndef HTSISAM_NO_DROPIN

/* sam_itr_next / sam_itr_multi_next are static inline functions in modern
 * htslib (1.10+). #define-shadowing redirects all callers in this TU. */
#undef sam_itr_next
#undef sam_itr_multi_next
#undef bam_itr_next
#undef bam_itr_querys
#undef bam_itr_queryi

#define hts_open            htsisam_hts_open
#define hts_open_format     htsisam_hts_open_format
#define sam_open            htsisam_sam_open
#define sam_close           htsisam_sam_close
#define hts_close           htsisam_hts_close

#define sam_hdr_read        htsisam_sam_hdr_read
#define sam_hdr_write       htsisam_sam_hdr_write

#define sam_read1           htsisam_sam_read1
#define sam_write1          htsisam_sam_write1

#define sam_index_load      htsisam_sam_index_load
#define sam_index_load2     htsisam_sam_index_load2
#define sam_index_load3     htsisam_sam_index_load3
#define hts_idx_destroy     htsisam_hts_idx_destroy

#define sam_itr_querys      htsisam_sam_itr_querys
#define sam_itr_queryi      htsisam_sam_itr_queryi
#define sam_itr_next        htsisam_sam_itr_next
#define sam_itr_multi_next  htsisam_sam_itr_next
#define hts_itr_destroy     htsisam_hts_itr_destroy

#define bam_itr_querys      htsisam_sam_itr_querys
#define bam_itr_queryi      htsisam_sam_itr_queryi
#define bam_itr_next        htsisam_sam_itr_next

#endif /* !HTSISAM_NO_DROPIN */

#endif /* HTSISAM_SAM_H */
