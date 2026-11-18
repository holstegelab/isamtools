#ifndef ISAM_STORE_H
#define ISAM_STORE_H

#include "common.h"
#include "backend_lmdb.h"

typedef struct cram_slice_cache cram_slice_cache_t;
typedef struct cram_decoder_pool cram_decoder_pool_t;

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
    cram_decoder_pool_t *decoders;
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
    lmdb_tag_fingerprint_t *tag_fingerprints;
    size_t   n_tag_fingerprints;
    int      has_tag_fingerprints;
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
int      store_record_tag_fingerprint(const store_record_t *rec,
                                      const char tag[2],
                                      uint64_t *fingerprint);

/* Fetch canonical payload (full BAM record) for the given store record.
 * Returns 0 on success, ISAM_ERR_* on failure. */
int      store_fetch_payload(store_t *st, const char *qname,
                              const store_record_t *rec, bam1_t *out);

/* Read-only payload, owned by the decoded CRAM cache or caller's scratch.
 * Valid until the next payload fetch on this store, or store_close. Prefetch
 * submission does not evict borrowed payloads. Never
 * mutate it. Copying store_fetch_payload remains available to API consumers. */
int      store_borrow_payload(store_t *st, const store_record_t *rec,
                               bam1_t *scratch, const bam1_t **payload);
/* Independent asynchronous container decoders. Prefetch never changes cache
 * residency or output order. Returns the examined prefix length; a shorter
 * prefix means the bounded queue is full and should be retried after fetching
 * the next demanded payload. The caller consumes records in coordinate order. */
void     store_set_decode_threads(store_t *st, int threads, const char *reference);
int      store_payload_cached(store_t *st, const store_record_t *rec);
size_t   store_prefetch_payloads(store_t *st, const store_record_t *records, size_t count);

/* Optional exact request plan for one coordinate traversal. Build before any
 * fetch/prefetch. Only requested records are retained from decoded containers;
 * release after each emitted/filtered record, when its borrowed payload is no
 * longer in use. The ordinary copy-based API does not enable this policy. */
void     store_payload_plan_begin(store_t *st);
void     store_payload_plan_add(store_t *st, const store_record_t *rec);
void     store_payload_plan_finish(store_t *st);
void     store_payload_release(store_t *st, const store_record_t *rec);

payload_cache_t *payload_cache_create(size_t max_bytes);
void             payload_cache_destroy(payload_cache_t *cache);
int              payload_cache_fetch(payload_cache_t *cache, store_t *st,
                                     const store_record_t *rec, bam1_t *out);

#endif
