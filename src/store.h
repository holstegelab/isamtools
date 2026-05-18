#ifndef ISAM_STORE_H
#define ISAM_STORE_H

#include "common.h"
#include "backend_lmdb.h"

typedef struct store {
    char       *store_path;
    char       *index_path;
    lmdb_backend_t *lmdb;
    samFile    *sam;           /* canonical BAM/CRAM */
    sam_hdr_t  *hdr;
    int         is_bam;        /* 1 for BAM, 0 for CRAM (or other) */
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
