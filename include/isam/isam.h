/*
 * libisam: compact, indexed, samtools-like alignment views backed by a
 * canonical BAM/CRAM read store. See DESIGN.md.
 *
 * v1 surface: enough to open an alignment .isam, read its header, iterate
 * records (sequential or by region) as fully reconstructed bam1_t records.
 */
#ifndef ISAM_H
#define ISAM_H

#include <stddef.h>
#include <stdint.h>
#include <htslib/sam.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#  define ISAM_EXPORT __declspec(dllexport)
#else
#  define ISAM_EXPORT __attribute__((visibility("default")))
#endif

#define ISAM_OK                  0
#define ISAM_EOF                -1
#define ISAM_ERR_IO             -2
#define ISAM_ERR_BACKEND        -3
#define ISAM_ERR_FORMAT         -4
#define ISAM_ERR_MISSING_STORE  -5
#define ISAM_ERR_MISSING_QNAME  -6
#define ISAM_ERR_STALE_INDEX    -7
#define ISAM_ERR_NO_INDEX       -8
#define ISAM_ERR_UNSUPPORTED    -9

typedef struct isamFile  isamFile;
typedef struct isamIndex isamIndex;
typedef struct isamItr   isamItr;

typedef struct isamOptions {
    const char *store_path_override;
    const char *store_index_path_override;
    const char *reference_fasta;
    int  load_store_annotations;
    int  load_view_annotations;
    int  tag_conflict_policy;     /* 0=view-wins, 1=store-wins, 2=error, 3=warn */
    size_t payload_cache_bytes;
    int  n_threads;
} isamOptions;

ISAM_EXPORT isamFile  *isam_open(const char *path, const char *mode);
ISAM_EXPORT isamFile  *isam_open_ext(const char *path, const char *mode,
                                     const isamOptions *opts);
ISAM_EXPORT int        isam_close(isamFile *fp);

ISAM_EXPORT sam_hdr_t *isam_hdr_read(isamFile *fp);
ISAM_EXPORT int        isam_hdr_write(isamFile *fp, const sam_hdr_t *hdr);

ISAM_EXPORT int        isam_read1(isamFile *fp, sam_hdr_t *hdr, bam1_t *b);

ISAM_EXPORT isamIndex *isam_index_load(const char *path);
ISAM_EXPORT void       isam_index_destroy(isamIndex *idx);

ISAM_EXPORT isamItr   *isam_itr_querys(isamIndex *idx, sam_hdr_t *hdr, const char *region);
ISAM_EXPORT isamItr   *isam_itr_queryi(isamIndex *idx, int tid,
                                       hts_pos_t beg, hts_pos_t end);
ISAM_EXPORT isamItr   *isam_itr_query_name(isamIndex *idx, const char *qname);
ISAM_EXPORT isamItr   *isam_itr_query_subset(isamIndex *idx, const char *subset_name);
ISAM_EXPORT isamItr   *isam_itr_query_region_subset(isamIndex *idx,
                                                    sam_hdr_t *hdr,
                                                    const char *region,
                                                    const char *subset_name);
ISAM_EXPORT int        isam_itr_next(isamFile *fp, isamItr *itr, bam1_t *b);
ISAM_EXPORT void       isam_itr_destroy(isamItr *itr);

ISAM_EXPORT const char *isam_strerror(int errcode);
ISAM_EXPORT int         isam_errno(isamFile *fp);

#ifdef __cplusplus
}
#endif

#endif /* ISAM_H */
