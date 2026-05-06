#ifndef ISAM_BACKEND_LMDB_H
#define ISAM_BACKEND_LMDB_H

#ifdef ISAM_ENABLE_LMDB

#include <stdint.h>
#include <stddef.h>
#include <lmdb.h>

typedef struct lmdb_backend {
    MDB_env *env;
    MDB_txn *write_txn;
    MDB_dbi meta;
    MDB_dbi canonical_by_id;
    MDB_dbi canonical_by_qname;
    MDB_dbi read_by_id;
    MDB_dbi read_by_qname;
    MDB_dbi aln_by_id;
    MDB_dbi aln_import_by_order;
    MDB_dbi aln_by_coord;
    MDB_dbi aln_unmapped;
    MDB_dbi aln_by_read;
    MDB_dbi aln_by_bin;
    char *path;
    int readonly;
} lmdb_backend_t;

typedef struct lmdb_store_record {
    int64_t read_id;
    char *qname;
    int64_t qname_hash;
    int tid;
    int64_t pos0;
    int64_t end0;
    int flag;
    int is_unmapped;
    int64_t bam_voffset;
    int has_voffset;
} lmdb_store_record_t;

typedef struct lmdb_view_record {
    int64_t aln_id;
    int64_t store_read_id;
    int flag;
    int tid;
    int64_t pos0;
    int64_t end0;
    int mapq;
    int n_cigar;
    uint32_t *cigar;
    int mate_tid;
    int64_t mate_pos0;
    int64_t tlen;
    uint8_t *aux;
    int aux_len;
    int64_t record_order;
    int is_unmapped;
} lmdb_view_record_t;

typedef enum lmdb_view_iter_kind {
    LMDB_VIEW_ITER_IMPORT = 0,
    LMDB_VIEW_ITER_COORD = 1,
    LMDB_VIEW_ITER_NAME = 2,
    LMDB_VIEW_ITER_REGION = 3
} lmdb_view_iter_kind_t;

typedef struct lmdb_view_iter lmdb_view_iter_t;

lmdb_backend_t *lmdb_backend_open(const char *path, int readonly,
                                  size_t map_size_bytes);
void            lmdb_backend_close(lmdb_backend_t *b);
void            lmdb_backend_begin(lmdb_backend_t *b);
void            lmdb_backend_commit(lmdb_backend_t *b);
void            lmdb_backend_abort(lmdb_backend_t *b);
void            lmdb_backend_meta_set(lmdb_backend_t *b,
                                      const char *key,
                                      const char *value);
char           *lmdb_backend_meta_get(lmdb_backend_t *b, const char *key);
int             lmdb_backend_meta_get_i64(lmdb_backend_t *b,
                                          const char *key,
                                          int64_t *out);

void            lmdb_store_put_meta(lmdb_backend_t *b,
                                    const char *key,
                                    const char *value);
void            lmdb_store_put_meta_i64(lmdb_backend_t *b,
                                        const char *key,
                                        int64_t value);
void            lmdb_store_put_record(lmdb_backend_t *b,
                                      const lmdb_store_record_t *rec);
int             lmdb_store_lookup_qname(lmdb_backend_t *b,
                                        int64_t qname_hash,
                                        const char *qname,
                                        lmdb_store_record_t *out);
int             lmdb_store_lookup_read_id(lmdb_backend_t *b,
                                          int64_t read_id,
                                          lmdb_store_record_t *out);
void            lmdb_store_record_clear(lmdb_store_record_t *rec);

void            lmdb_view_put_meta(lmdb_backend_t *b,
                                   const char *key,
                                   const char *value);
void            lmdb_view_put_meta_i64(lmdb_backend_t *b,
                                       const char *key,
                                       int64_t value);
void            lmdb_view_delete_meta(lmdb_backend_t *b, const char *key);
void            lmdb_view_put_read(lmdb_backend_t *b,
                                   int64_t store_read_id,
                                   const char *qname,
                                   int64_t qname_hash);
int             lmdb_view_lookup_qname(lmdb_backend_t *b,
                                       int64_t qname_hash,
                                       const char *qname,
                                       int64_t *store_read_id);
void            lmdb_view_put_import_record(lmdb_backend_t *b,
                                            const lmdb_view_record_t *rec);
void            lmdb_view_clear_sorted(lmdb_backend_t *b,
                                       int clear_coord,
                                       int clear_interval);
void            lmdb_view_index_record(lmdb_backend_t *b,
                                       const lmdb_view_record_t *rec,
                                       int index_coord,
                                       int index_interval);
lmdb_view_iter_t *lmdb_view_iter_open(lmdb_backend_t *b,
                                      lmdb_view_iter_kind_t kind,
                                      int tid,
                                      int64_t beg,
                                      int64_t end);
int             lmdb_view_iter_next(lmdb_view_iter_t *it,
                                    lmdb_view_record_t *out);
void            lmdb_view_iter_close(lmdb_view_iter_t *it);
void            lmdb_view_record_clear(lmdb_view_record_t *rec);

#endif

#endif
