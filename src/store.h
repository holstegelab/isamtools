#ifndef ISAM_STORE_H
#define ISAM_STORE_H

#include "common.h"
#include "backend_lmdb.h"

typedef struct cram_slice_cache cram_slice_cache_t;

typedef struct store {
    char       *store_path;
    char       *index_path;
    lmdb_backend_t *lmdb;
    samFile    *sam;           /* canonical BAM/CRAM */
    sam_hdr_t  *hdr;
    int         is_bam;        /* 1 for BAM, 0 for CRAM (or other) */
    /* CRAM canonical: container table (offset -> record count) + decoded-slice
     * LRU cache, so payload fetch is by container locator, decoding each
     * container's slice at most once while it stays resident. */
    int64_t    *cram_coff;     /* container byte offsets, ascending; NULL for BAM */
    int32_t    *cram_nrec;     /* records per container */
    int64_t     cram_nc;       /* number of containers */
    cram_slice_cache_t *slice_cache;
} store_t;

typedef struct store_record {
    int64_t read_id;
    char   *qname;             /* owned by caller; free with store_record_clear */
    int     tid;
    int64_t pos0;
    int64_t end0;
    int     flag;
    int     is_unmapped;
    int64_t bam_voffset;       /* may be -1 if absent */
    int     has_voffset;
    int64_t cram_coffset;      /* CRAM container byte offset, -1 if absent */
    int32_t cram_within;       /* record index within the container */
    int     has_cram_loc;
} store_record_t;

typedef struct payload_cache payload_cache_t;

/* Open a canonical store + its .ssam index (read-only). */
store_t *store_open(const char *store_path, const char *index_path);
void     store_close(store_t *st);

/* Lookup QNAME in the store index. Returns 1 on hit, 0 on miss. */
int      store_has_name_index(store_t *st);
int      store_lookup_qname(store_t *st, const char *qname, store_record_t *out);
int      store_lookup_read_id(store_t *st, int64_t read_id, store_record_t *out);
void     store_record_clear(store_record_t *rec);

/* Fetch canonical payload (full BAM record) for the given store record.
 * Returns 0 on success, ISAM_ERR_* on failure. */
int      store_fetch_payload(store_t *st, const char *qname,
                              const store_record_t *rec, bam1_t *out);

payload_cache_t *payload_cache_create(size_t max_bytes);
void             payload_cache_destroy(payload_cache_t *cache);
int              payload_cache_fetch(payload_cache_t *cache, store_t *st,
                                     const store_record_t *rec, bam1_t *out);

#endif
