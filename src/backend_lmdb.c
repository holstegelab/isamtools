#include "common.h"
#include "backend_lmdb.h"

#ifdef ISAM_ENABLE_LMDB

#include <sys/stat.h>
#include <errno.h>

static void lmdb_die(int rc, const char *op, const char *path) {
    die("lmdb %s failed for '%s': %s", op, path, mdb_strerror(rc));
}

static void put_u16(unsigned char **p, uint16_t v) {
    *(*p)++ = (unsigned char)(v >> 8);
    *(*p)++ = (unsigned char)v;
}

static void put_u32(unsigned char **p, uint32_t v) {
    for (int i = 3; i >= 0; --i) *(*p)++ = (unsigned char)(v >> (i * 8));
}

static void put_u64(unsigned char **p, uint64_t v) {
    for (int i = 7; i >= 0; --i) *(*p)++ = (unsigned char)(v >> (i * 8));
}

static uint16_t get_u16(const unsigned char **p) {
    uint16_t v = ((uint16_t)(*p)[0] << 8) | (*p)[1];
    *p += 2;
    return v;
}

static uint32_t get_u32(const unsigned char **p) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v = (v << 8) | *(*p)++;
    return v;
}

static uint64_t get_u64(const unsigned char **p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 8) | *(*p)++;
    return v;
}

static void be64_key(int64_t value, unsigned char key[8]) {
    unsigned char *p = key;
    put_u64(&p, (uint64_t)value);
}

static int64_t decode_i64_value(const MDB_val *val, const char *what) {
    if (val->mv_size != 8) die("corrupt lmdb %s value", what);
    const unsigned char *p = val->mv_data;
    return (int64_t)get_u64(&p);
}

static MDB_val qname_key(int64_t hash, const char *qname) {
    size_t n = strlen(qname);
    unsigned char *buf = xmalloc(8 + n);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)hash);
    memcpy(p, qname, n);
    MDB_val key = { 8 + n, buf };
    return key;
}

static MDB_txn *write_txn(lmdb_backend_t *b, int *owned) {
    if (b->readonly) die("lmdb write on readonly backend '%s'", b->path);
    if (b->write_txn) {
        *owned = 0;
        return b->write_txn;
    }
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, 0, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    *owned = 1;
    return txn;
}

static void write_txn_finish(lmdb_backend_t *b, MDB_txn *txn, int owned) {
    if (!owned) return;
    int rc = mdb_txn_commit(txn);
    if (rc) lmdb_die(rc, "txn_commit", b->path);
}

static void dbi_open_required(MDB_txn *txn, const char *path,
                              const char *name, unsigned int flags,
                              MDB_dbi *out) {
    int rc = mdb_dbi_open(txn, name, flags, out);
    if (rc) lmdb_die(rc, name, path);
}

static void dbi_open_optional(MDB_txn *txn, const char *path,
                              const char *name, unsigned int flags,
                              MDB_dbi *out) {
    int rc = mdb_dbi_open(txn, name, flags, out);
    if (rc == MDB_NOTFOUND) {
        *out = 0;
        return;
    }
    if (rc) lmdb_die(rc, name, path);
}

static MDB_val order_key(int64_t record_order, int64_t aln_id) {
    unsigned char *buf = xmalloc(16);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)record_order);
    put_u64(&p, (uint64_t)aln_id);
    MDB_val key = { 16, buf };
    return key;
}

static MDB_val coord_key(int tid, int64_t pos0,
                         int64_t record_order, int64_t aln_id) {
    unsigned char *buf = xmalloc(32);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)tid);
    put_u64(&p, (uint64_t)pos0);
    put_u64(&p, (uint64_t)record_order);
    put_u64(&p, (uint64_t)aln_id);
    MDB_val key = { 32, buf };
    return key;
}

static MDB_val read_aln_key(int64_t store_read_id,
                            int64_t record_order,
                            int64_t aln_id) {
    unsigned char *buf = xmalloc(24);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)store_read_id);
    put_u64(&p, (uint64_t)record_order);
    put_u64(&p, (uint64_t)aln_id);
    MDB_val key = { 24, buf };
    return key;
}

static MDB_val bin_key(int tid, int64_t bin,
                       int64_t pos0, int64_t record_order, int64_t aln_id) {
    unsigned char *buf = xmalloc(40);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)tid);
    put_u64(&p, (uint64_t)bin);
    put_u64(&p, (uint64_t)pos0);
    put_u64(&p, (uint64_t)record_order);
    put_u64(&p, (uint64_t)aln_id);
    MDB_val key = { 40, buf };
    return key;
}

static int key_has_prefix(const MDB_val *key, const unsigned char *prefix,
                          size_t prefix_len) {
    return key->mv_size >= prefix_len &&
           memcmp(key->mv_data, prefix, prefix_len) == 0;
}

static MDB_val encode_record(const lmdb_store_record_t *rec) {
    size_t qn = strlen(rec->qname);
    if (qn > UINT16_MAX) die("qname too long for lmdb record");
    size_t n = 1 + 2 + qn + 8 + 8 + 4 + 8 + 8 + 4 + 4 + 4 + 8 + 4;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    put_u16(&p, (uint16_t)qn);
    memcpy(p, rec->qname, qn); p += qn;
    put_u64(&p, (uint64_t)rec->read_id);
    put_u64(&p, (uint64_t)rec->qname_hash);
    put_u32(&p, (uint32_t)rec->tid);
    put_u64(&p, (uint64_t)rec->pos0);
    put_u64(&p, (uint64_t)rec->end0);
    put_u32(&p, (uint32_t)rec->flag);
    put_u32(&p, (uint32_t)rec->is_unmapped);
    put_u32(&p, (uint32_t)rec->has_voffset);
    put_u64(&p, (uint64_t)rec->bam_voffset);
    put_u32(&p, 0);
    MDB_val val = { n, buf };
    return val;
}

static void decode_record(MDB_val *val, lmdb_store_record_t *out) {
    memset(out, 0, sizeof *out);
    const unsigned char *p = val->mv_data;
    const unsigned char *end = p + val->mv_size;
    if (p >= end || *p++ != 1) die("unsupported lmdb canonical record version");
    uint16_t qn = get_u16(&p);
    if (p + qn > end) die("truncated lmdb canonical record");
    out->qname = xmalloc((size_t)qn + 1);
    memcpy(out->qname, p, qn);
    out->qname[qn] = 0;
    p += qn;
    out->read_id = (int64_t)get_u64(&p);
    out->qname_hash = (int64_t)get_u64(&p);
    out->tid = (int32_t)get_u32(&p);
    out->pos0 = (int64_t)get_u64(&p);
    out->end0 = (int64_t)get_u64(&p);
    out->flag = (int32_t)get_u32(&p);
    out->is_unmapped = (int32_t)get_u32(&p);
    out->has_voffset = (int32_t)get_u32(&p);
    out->bam_voffset = (int64_t)get_u64(&p);
    (void)get_u32(&p);
}

static MDB_val encode_view_read(int64_t store_read_id,
                                const char *qname,
                                int64_t hash) {
    size_t qn = strlen(qname);
    if (qn > UINT16_MAX) die("qname too long for lmdb view read");
    size_t n = 1 + 2 + qn + 8 + 8;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    put_u16(&p, (uint16_t)qn);
    memcpy(p, qname, qn); p += qn;
    put_u64(&p, (uint64_t)store_read_id);
    put_u64(&p, (uint64_t)hash);
    MDB_val val = { n, buf };
    return val;
}

static MDB_val encode_view_record(const lmdb_view_record_t *rec) {
    if (rec->n_cigar < 0 || rec->aux_len < 0)
        die("negative blob length in lmdb view record");
    size_t cigar_bytes = (size_t)rec->n_cigar * 4;
    size_t aux_bytes = (size_t)rec->aux_len;
    size_t n = 1 + 8 + 8 + 4 + 4 + 8 + 8 + 4 + 4 + 4 + 8 + 8 +
               8 + 4 + 4 + 8 + 4 + cigar_bytes + aux_bytes;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    put_u64(&p, (uint64_t)rec->aln_id);
    put_u64(&p, (uint64_t)rec->store_read_id);
    put_u32(&p, (uint32_t)rec->flag);
    put_u32(&p, (uint32_t)rec->tid);
    put_u64(&p, (uint64_t)rec->pos0);
    put_u64(&p, (uint64_t)rec->end0);
    put_u32(&p, (uint32_t)rec->mapq);
    put_u32(&p, (uint32_t)rec->n_cigar);
    put_u32(&p, (uint32_t)rec->mate_tid);
    put_u64(&p, (uint64_t)rec->mate_pos0);
    put_u64(&p, (uint64_t)rec->tlen);
    put_u64(&p, (uint64_t)rec->record_order);
    put_u32(&p, (uint32_t)rec->is_unmapped);
    put_u32(&p, (uint32_t)rec->aux_len);
    put_u64(&p, 0);
    put_u32(&p, 0);
    if (cigar_bytes) { memcpy(p, rec->cigar, cigar_bytes); p += cigar_bytes; }
    if (aux_bytes) { memcpy(p, rec->aux, aux_bytes); p += aux_bytes; }
    MDB_val val = { n, buf };
    return val;
}

static void decode_view_record(const MDB_val *val, lmdb_view_record_t *out) {
    memset(out, 0, sizeof *out);
    const unsigned char *p = val->mv_data;
    const unsigned char *end = p + val->mv_size;
    if (p >= end || *p++ != 1) die("unsupported lmdb view record version");
    out->aln_id = (int64_t)get_u64(&p);
    out->store_read_id = (int64_t)get_u64(&p);
    out->flag = (int32_t)get_u32(&p);
    out->tid = (int32_t)get_u32(&p);
    out->pos0 = (int64_t)get_u64(&p);
    out->end0 = (int64_t)get_u64(&p);
    out->mapq = (int32_t)get_u32(&p);
    out->n_cigar = (int32_t)get_u32(&p);
    out->mate_tid = (int32_t)get_u32(&p);
    out->mate_pos0 = (int64_t)get_u64(&p);
    out->tlen = (int64_t)get_u64(&p);
    out->record_order = (int64_t)get_u64(&p);
    out->is_unmapped = (int32_t)get_u32(&p);
    out->aux_len = (int32_t)get_u32(&p);
    (void)get_u64(&p);
    (void)get_u32(&p);
    if (out->n_cigar < 0 || out->aux_len < 0)
        die("corrupt lmdb view record blob length");
    size_t cigar_bytes = (size_t)out->n_cigar * 4;
    size_t aux_bytes = (size_t)out->aux_len;
    if (p + cigar_bytes + aux_bytes > end) die("truncated lmdb view record");
    if (cigar_bytes) {
        out->cigar = xmalloc(cigar_bytes);
        memcpy(out->cigar, p, cigar_bytes);
        p += cigar_bytes;
    }
    if (aux_bytes) {
        out->aux = xmalloc(aux_bytes);
        memcpy(out->aux, p, aux_bytes);
    }
}

static int lmdb_view_get_record_txn(lmdb_backend_t *b,
                                    MDB_txn *txn,
                                    int64_t aln_id,
                                    lmdb_view_record_t *out) {
    unsigned char idbuf[8];
    be64_key(aln_id, idbuf);
    MDB_val key = { sizeof idbuf, idbuf };
    MDB_val val = { 0, NULL };
    int rc = mdb_get(txn, b->aln_by_id, &key, &val);
    if (rc == MDB_NOTFOUND) return 0;
    if (rc) lmdb_die(rc, "get(aln_by_id)", b->path);
    decode_view_record(&val, out);
    return 1;
}

lmdb_backend_t *lmdb_backend_open(const char *path, int readonly,
                                  size_t map_size_bytes) {
    lmdb_backend_t *b = xcalloc(1, sizeof *b);
    b->path = xstrdup(path);
    b->readonly = readonly;

    int rc = mdb_env_create(&b->env);
    if (rc) lmdb_die(rc, "env_create", path);

    rc = mdb_env_set_maxdbs(b->env, 64);
    if (rc) lmdb_die(rc, "env_set_maxdbs", path);

    rc = mdb_env_set_mapsize(b->env,
                             map_size_bytes ? map_size_bytes
                                            : (size_t)64 * 1024 * 1024 * 1024);
    if (rc) lmdb_die(rc, "env_set_mapsize", path);

    int flags = MDB_NOSUBDIR | (readonly ? MDB_RDONLY : 0);
    rc = mdb_env_open(b->env, path, flags, 0664);
    if (rc) lmdb_die(rc, "env_open", path);

    MDB_txn *txn = NULL;
    rc = mdb_txn_begin(b->env, NULL, readonly ? MDB_RDONLY : 0, &txn);
    if (rc) lmdb_die(rc, "txn_begin", path);

    dbi_open_required(txn, path, "meta", readonly ? 0 : MDB_CREATE, &b->meta);

    unsigned int dbi_flags = readonly ? 0 : MDB_CREATE;
    dbi_open_optional(txn, path, "canonical_by_id", dbi_flags,
                      &b->canonical_by_id);
    dbi_open_optional(txn, path, "canonical_by_qname", dbi_flags,
                      &b->canonical_by_qname);
    dbi_open_optional(txn, path, "read_by_id", dbi_flags, &b->read_by_id);
    dbi_open_optional(txn, path, "read_by_qname", dbi_flags, &b->read_by_qname);
    dbi_open_optional(txn, path, "aln_by_id", dbi_flags, &b->aln_by_id);
    dbi_open_optional(txn, path, "aln_import_by_order", dbi_flags,
                      &b->aln_import_by_order);
    dbi_open_optional(txn, path, "aln_by_coord", dbi_flags, &b->aln_by_coord);
    dbi_open_optional(txn, path, "aln_unmapped", dbi_flags, &b->aln_unmapped);
    dbi_open_optional(txn, path, "aln_by_read", dbi_flags, &b->aln_by_read);
    dbi_open_optional(txn, path, "aln_by_bin", dbi_flags, &b->aln_by_bin);

    rc = mdb_txn_commit(txn);
    if (rc) lmdb_die(rc, "txn_commit", path);
    return b;
}

void lmdb_backend_close(lmdb_backend_t *b) {
    if (!b) return;
    if (b->write_txn) lmdb_backend_abort(b);
    if (b->env) {
        mdb_dbi_close(b->env, b->meta);
        if (b->canonical_by_id) mdb_dbi_close(b->env, b->canonical_by_id);
        if (b->canonical_by_qname) mdb_dbi_close(b->env, b->canonical_by_qname);
        if (b->read_by_id) mdb_dbi_close(b->env, b->read_by_id);
        if (b->read_by_qname) mdb_dbi_close(b->env, b->read_by_qname);
        if (b->aln_by_id) mdb_dbi_close(b->env, b->aln_by_id);
        if (b->aln_import_by_order) mdb_dbi_close(b->env, b->aln_import_by_order);
        if (b->aln_by_coord) mdb_dbi_close(b->env, b->aln_by_coord);
        if (b->aln_unmapped) mdb_dbi_close(b->env, b->aln_unmapped);
        if (b->aln_by_read) mdb_dbi_close(b->env, b->aln_by_read);
        if (b->aln_by_bin) mdb_dbi_close(b->env, b->aln_by_bin);
        mdb_env_close(b->env);
    }
    free(b->path);
    free(b);
}

void lmdb_backend_begin(lmdb_backend_t *b) {
    if (b->readonly) die("lmdb begin on readonly backend '%s'", b->path);
    if (b->write_txn) die("nested lmdb write transaction on '%s'", b->path);
    int rc = mdb_txn_begin(b->env, NULL, 0, &b->write_txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
}

void lmdb_backend_commit(lmdb_backend_t *b) {
    if (!b->write_txn) return;
    MDB_txn *txn = b->write_txn;
    b->write_txn = NULL;
    int rc = mdb_txn_commit(txn);
    if (rc) lmdb_die(rc, "txn_commit", b->path);
}

void lmdb_backend_abort(lmdb_backend_t *b) {
    if (!b->write_txn) return;
    MDB_txn *txn = b->write_txn;
    b->write_txn = NULL;
    mdb_txn_abort(txn);
}

void lmdb_backend_meta_set(lmdb_backend_t *b,
                           const char *key,
                           const char *value) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);

    MDB_val k = { strlen(key), (void *)key };
    MDB_val v = { strlen(value), (void *)value };
    int rc = mdb_put(txn, b->meta, &k, &v, 0);
    if (rc) lmdb_die(rc, "put(meta)", b->path);

    write_txn_finish(b, txn, owned);
}

char *lmdb_backend_meta_get(lmdb_backend_t *b, const char *key) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);

    MDB_val k = { strlen(key), (void *)key };
    MDB_val v = { 0, NULL };
    rc = mdb_get(txn, b->meta, &k, &v);
    if (rc == MDB_NOTFOUND) {
        mdb_txn_abort(txn);
        return NULL;
    }
    if (rc) lmdb_die(rc, "get(meta)", b->path);

    char *out = xmalloc(v.mv_size + 1);
    memcpy(out, v.mv_data, v.mv_size);
    out[v.mv_size] = 0;
    mdb_txn_abort(txn);
    return out;
}

int lmdb_backend_meta_get_i64(lmdb_backend_t *b, const char *key, int64_t *out) {
    char *s = lmdb_backend_meta_get(b, key);
    if (!s) return 0;
    *out = (int64_t)strtoll(s, NULL, 10);
    free(s);
    return 1;
}

void lmdb_store_put_meta(lmdb_backend_t *b,
                         const char *key,
                         const char *value) {
    lmdb_backend_meta_set(b, key, value);
}

void lmdb_store_put_meta_i64(lmdb_backend_t *b,
                             const char *key,
                             int64_t value) {
    char buf[32];
    snprintf(buf, sizeof buf, "%lld", (long long)value);
    lmdb_backend_meta_set(b, key, buf);
}

void lmdb_store_put_record(lmdb_backend_t *b,
                           const lmdb_store_record_t *rec) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);

    unsigned char idbuf[8];
    be64_key(rec->read_id, idbuf);
    MDB_val idkey = { sizeof idbuf, idbuf };
    MDB_val val = encode_record(rec);
    int rc = mdb_put(txn, b->canonical_by_id, &idkey, &val, MDB_NOOVERWRITE);
    if (rc) lmdb_die(rc, "put(canonical_by_id)", b->path);

    MDB_val qkey = qname_key(rec->qname_hash, rec->qname);
    unsigned char ridbuf[8];
    unsigned char *p = ridbuf;
    put_u64(&p, (uint64_t)rec->read_id);
    MDB_val ridval = { sizeof ridbuf, ridbuf };
    rc = mdb_put(txn, b->canonical_by_qname, &qkey, &ridval, MDB_NOOVERWRITE);
    if (rc) lmdb_die(rc, "put(canonical_by_qname)", b->path);

    free(val.mv_data);
    free(qkey.mv_data);
    write_txn_finish(b, txn, owned);
}

int lmdb_store_lookup_read_id(lmdb_backend_t *b,
                              int64_t read_id,
                              lmdb_store_record_t *out) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);

    unsigned char idbuf[8];
    be64_key(read_id, idbuf);
    MDB_val key = { sizeof idbuf, idbuf };
    MDB_val val = { 0, NULL };
    rc = mdb_get(txn, b->canonical_by_id, &key, &val);
    if (rc == MDB_NOTFOUND) {
        mdb_txn_abort(txn);
        return 0;
    }
    if (rc) lmdb_die(rc, "get(canonical_by_id)", b->path);
    decode_record(&val, out);
    mdb_txn_abort(txn);
    return 1;
}

int lmdb_store_lookup_qname(lmdb_backend_t *b,
                            int64_t hash,
                            const char *qname,
                            lmdb_store_record_t *out) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);

    MDB_val key = qname_key(hash, qname);
    MDB_val val = { 0, NULL };
    rc = mdb_get(txn, b->canonical_by_qname, &key, &val);
    free(key.mv_data);
    if (rc == MDB_NOTFOUND) {
        mdb_txn_abort(txn);
        return 0;
    }
    if (rc) lmdb_die(rc, "get(canonical_by_qname)", b->path);
    const unsigned char *p = val.mv_data;
    int64_t read_id = (int64_t)get_u64(&p);
    mdb_txn_abort(txn);
    return lmdb_store_lookup_read_id(b, read_id, out);
}

void lmdb_store_record_clear(lmdb_store_record_t *rec) {
    if (!rec) return;
    free(rec->qname);
    memset(rec, 0, sizeof *rec);
}

void lmdb_view_put_meta(lmdb_backend_t *b,
                        const char *key,
                        const char *value) {
    lmdb_backend_meta_set(b, key, value);
}

void lmdb_view_put_meta_i64(lmdb_backend_t *b,
                            const char *key,
                            int64_t value) {
    char buf[32];
    snprintf(buf, sizeof buf, "%lld", (long long)value);
    lmdb_backend_meta_set(b, key, buf);
}

void lmdb_view_delete_meta(lmdb_backend_t *b, const char *key) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    MDB_val k = { strlen(key), (void *)key };
    int rc = mdb_del(txn, b->meta, &k, NULL);
    if (rc && rc != MDB_NOTFOUND) lmdb_die(rc, "del(meta)", b->path);
    write_txn_finish(b, txn, owned);
}

void lmdb_view_put_read(lmdb_backend_t *b,
                        int64_t store_read_id,
                        const char *qname,
                        int64_t hash) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);

    unsigned char idbuf[8];
    be64_key(store_read_id, idbuf);
    MDB_val idkey = { sizeof idbuf, idbuf };
    MDB_val val = encode_view_read(store_read_id, qname, hash);
    int rc = mdb_put(txn, b->read_by_id, &idkey, &val, 0);
    if (rc) lmdb_die(rc, "put(read_by_id)", b->path);

    MDB_val qkey = qname_key(hash, qname);
    unsigned char ridbuf[8];
    unsigned char *p = ridbuf;
    put_u64(&p, (uint64_t)store_read_id);
    MDB_val ridval = { sizeof ridbuf, ridbuf };
    rc = mdb_put(txn, b->read_by_qname, &qkey, &ridval, 0);
    if (rc) lmdb_die(rc, "put(read_by_qname)", b->path);

    free(val.mv_data);
    free(qkey.mv_data);
    write_txn_finish(b, txn, owned);
}

int lmdb_view_lookup_qname(lmdb_backend_t *b,
                           int64_t hash,
                           const char *qname,
                           int64_t *store_read_id) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    MDB_val key = qname_key(hash, qname);
    MDB_val val = { 0, NULL };
    rc = mdb_get(txn, b->read_by_qname, &key, &val);
    free(key.mv_data);
    if (rc == MDB_NOTFOUND) {
        mdb_txn_abort(txn);
        return 0;
    }
    if (rc) lmdb_die(rc, "get(read_by_qname)", b->path);
    const unsigned char *p = val.mv_data;
    *store_read_id = (int64_t)get_u64(&p);
    mdb_txn_abort(txn);
    return 1;
}

void lmdb_view_put_import_record(lmdb_backend_t *b,
                                 const lmdb_view_record_t *rec) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    MDB_val val = encode_view_record(rec);

    unsigned char idbuf[8];
    be64_key(rec->aln_id, idbuf);
    MDB_val idkey = { sizeof idbuf, idbuf };
    int rc = mdb_put(txn, b->aln_by_id, &idkey, &val, MDB_NOOVERWRITE);
    if (rc) lmdb_die(rc, "put(aln_by_id)", b->path);

    unsigned char refbuf[8];
    unsigned char *p = refbuf;
    put_u64(&p, (uint64_t)rec->aln_id);
    MDB_val refval = { sizeof refbuf, refbuf };

    MDB_val okey = order_key(rec->record_order, rec->aln_id);
    rc = mdb_put(txn, b->aln_import_by_order, &okey, &refval, 0);
    if (rc) lmdb_die(rc, "put(aln_import_by_order)", b->path);

    MDB_val rkey = read_aln_key(rec->store_read_id, rec->record_order,
                                rec->aln_id);
    rc = mdb_put(txn, b->aln_by_read, &rkey, &refval, 0);
    if (rc) lmdb_die(rc, "put(aln_by_read)", b->path);

    free(okey.mv_data);
    free(rkey.mv_data);
    free(val.mv_data);
    write_txn_finish(b, txn, owned);
}

static void lmdb_drop_clear(lmdb_backend_t *b, MDB_txn *txn,
                            MDB_dbi dbi, const char *name) {
    int rc = mdb_drop(txn, dbi, 0);
    if (rc) lmdb_die(rc, name, b->path);
}

void lmdb_view_clear_sorted(lmdb_backend_t *b,
                            int clear_coord,
                            int clear_interval) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    if (clear_coord) {
        lmdb_drop_clear(b, txn, b->aln_by_coord, "drop(aln_by_coord)");
        lmdb_drop_clear(b, txn, b->aln_unmapped, "drop(aln_unmapped)");
    }
    if (clear_interval)
        lmdb_drop_clear(b, txn, b->aln_by_bin, "drop(aln_by_bin)");
    write_txn_finish(b, txn, owned);
}

void lmdb_view_index_record(lmdb_backend_t *b,
                            const lmdb_view_record_t *rec,
                            int index_coord,
                            int index_interval) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    unsigned char refbuf[8];
    unsigned char *p = refbuf;
    put_u64(&p, (uint64_t)rec->aln_id);
    MDB_val refval = { sizeof refbuf, refbuf };
    int mapped = !rec->is_unmapped && rec->tid >= 0 && rec->pos0 >= 0;

    if (index_coord) {
        MDB_val key = mapped
            ? coord_key(rec->tid, rec->pos0, rec->record_order, rec->aln_id)
            : order_key(rec->record_order, rec->aln_id);
        MDB_dbi dbi = mapped ? b->aln_by_coord : b->aln_unmapped;
        int rc = mdb_put(txn, dbi, &key, &refval, 0);
        if (rc) lmdb_die(rc, mapped ? "put(aln_by_coord)" : "put(aln_unmapped)",
                         b->path);
        free(key.mv_data);
    }

    if (index_interval && mapped) {
        const int64_t bin_size = 16384;
        int64_t end0 = rec->end0 > rec->pos0 ? rec->end0 : rec->pos0 + 1;
        int64_t first_bin = rec->pos0 / bin_size;
        int64_t last_bin = (end0 - 1) / bin_size;
        for (int64_t bin = first_bin; bin <= last_bin; ++bin) {
            MDB_val key = bin_key(rec->tid, bin, rec->pos0,
                                  rec->record_order, rec->aln_id);
            int rc = mdb_put(txn, b->aln_by_bin, &key, &refval, 0);
            if (rc) lmdb_die(rc, "put(aln_by_bin)", b->path);
            free(key.mv_data);
        }
    }

    write_txn_finish(b, txn, owned);
}

struct lmdb_view_iter {
    lmdb_backend_t *b;
    MDB_txn *txn;
    MDB_cursor *cursor;
    MDB_cursor *read_cursor;
    lmdb_view_iter_kind_t kind;
    int started;
    int phase;
    int tid;
    int64_t beg;
    int64_t end;
    int64_t bin;
    int64_t last_bin;
    unsigned char read_prefix[8];
    int have_read_prefix;
    int64_t *seen_aln_ids;
    size_t n_seen;
    size_t cap_seen;
};

static void seen_add(lmdb_view_iter_t *it, int64_t aln_id) {
    if (it->n_seen == it->cap_seen) {
        it->cap_seen = it->cap_seen ? it->cap_seen * 2 : 64;
        it->seen_aln_ids = xrealloc(it->seen_aln_ids,
                                    it->cap_seen * sizeof *it->seen_aln_ids);
    }
    it->seen_aln_ids[it->n_seen++] = aln_id;
}

static int seen_has(lmdb_view_iter_t *it, int64_t aln_id) {
    for (size_t i = 0; i < it->n_seen; ++i)
        if (it->seen_aln_ids[i] == aln_id) return 1;
    return 0;
}

lmdb_view_iter_t *lmdb_view_iter_open(lmdb_backend_t *b,
                                      lmdb_view_iter_kind_t kind,
                                      int tid,
                                      int64_t beg,
                                      int64_t end) {
    lmdb_view_iter_t *it = xcalloc(1, sizeof *it);
    it->b = b;
    it->kind = kind;
    it->tid = tid;
    it->beg = beg;
    it->end = end;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &it->txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);

    MDB_dbi dbi = b->aln_import_by_order;
    if (kind == LMDB_VIEW_ITER_COORD) dbi = b->aln_by_coord;
    else if (kind == LMDB_VIEW_ITER_REGION) {
        dbi = b->aln_by_bin;
        const int64_t bin_size = 16384;
        int64_t query_end = end > beg ? end : beg + 1;
        it->bin = beg / bin_size;
        it->last_bin = (query_end - 1) / bin_size;
    } else if (kind == LMDB_VIEW_ITER_NAME) {
        rc = mdb_cursor_open(it->txn, b->read_by_qname, &it->read_cursor);
        if (rc) lmdb_die(rc, "cursor(read_by_qname)", b->path);
        dbi = b->aln_by_read;
    }
    rc = mdb_cursor_open(it->txn, dbi, &it->cursor);
    if (rc) lmdb_die(rc, "cursor(view)", b->path);
    return it;
}

static int iter_cursor_next(lmdb_view_iter_t *it, MDB_val *key, MDB_val *val) {
    int rc = mdb_cursor_get(it->cursor, key, val,
                            it->started ? MDB_NEXT : MDB_FIRST);
    it->started = 1;
    if (rc == MDB_NOTFOUND) return 0;
    if (rc) lmdb_die(rc, "cursor_get", it->b->path);
    return 1;
}

static int iter_decode_ref(lmdb_view_iter_t *it,
                           const MDB_val *val,
                           lmdb_view_record_t *out,
                           const char *what) {
    int64_t aln_id = decode_i64_value(val, what);
    if (!lmdb_view_get_record_txn(it->b, it->txn, aln_id, out))
        die("lmdb alignment record '%lld' referenced by %s is missing in '%s'",
            (long long)aln_id, what, it->b->path);
    return 1;
}

static int iter_coord_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    if (iter_cursor_next(it, &key, &val))
        return iter_decode_ref(it, &val, out, "aln_by_coord");
    if (it->phase == 0) {
        mdb_cursor_close(it->cursor);
        it->cursor = NULL;
        int rc = mdb_cursor_open(it->txn, it->b->aln_unmapped, &it->cursor);
        if (rc) lmdb_die(rc, "cursor(aln_unmapped)", it->b->path);
        it->started = 0;
        it->phase = 1;
        if (iter_cursor_next(it, &key, &val))
            return iter_decode_ref(it, &val, out, "aln_unmapped");
    }
    return 0;
}

static int iter_name_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while (1) {
        if (it->have_read_prefix) {
            int rc = mdb_cursor_get(it->cursor, &key, &val,
                                    it->started ? MDB_NEXT : MDB_SET_RANGE);
            it->started = 1;
            if (rc && rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(aln_by_read)", it->b->path);
            if (!rc && key_has_prefix(&key, it->read_prefix, sizeof it->read_prefix))
                return iter_decode_ref(it, &val, out, "aln_by_read");
            it->have_read_prefix = 0;
            it->started = 0;
        }

        int rc = mdb_cursor_get(it->read_cursor, &key, &val,
                                it->phase ? MDB_NEXT : MDB_FIRST);
        it->phase = 1;
        if (rc == MDB_NOTFOUND) return 0;
        if (rc) lmdb_die(rc, "cursor(read_by_qname)", it->b->path);
        const unsigned char *p = val.mv_data;
        int64_t store_read_id = (int64_t)get_u64(&p);
        unsigned char *kp = it->read_prefix;
        put_u64(&kp, (uint64_t)store_read_id);
        key.mv_size = sizeof it->read_prefix;
        key.mv_data = it->read_prefix;
        it->have_read_prefix = 1;
        it->started = 0;
    }
}

static int iter_region_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    const int64_t bin_size = 16384;
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while (it->bin <= it->last_bin) {
        unsigned char prefix[16];
        unsigned char *p = prefix;
        put_u64(&p, (uint64_t)it->tid);
        put_u64(&p, (uint64_t)it->bin);

        int rc;
        if (!it->started) {
            key.mv_size = sizeof prefix;
            key.mv_data = prefix;
            rc = mdb_cursor_get(it->cursor, &key, &val, MDB_SET_RANGE);
            it->started = 1;
        } else {
            rc = mdb_cursor_get(it->cursor, &key, &val, MDB_NEXT);
        }
        if (rc == MDB_NOTFOUND) {
            ++it->bin;
            it->started = 0;
            continue;
        }
        if (rc) lmdb_die(rc, "cursor(aln_by_bin)", it->b->path);
        if (!key_has_prefix(&key, prefix, sizeof prefix)) {
            ++it->bin;
            it->started = 0;
            continue;
        }

        int64_t aln_id = decode_i64_value(&val, "aln_by_bin");
        lmdb_view_record_t rec;
        if (!lmdb_view_get_record_txn(it->b, it->txn, aln_id, &rec))
            die("lmdb alignment record '%lld' referenced by aln_by_bin is missing in '%s'",
                (long long)aln_id, it->b->path);
        int overlaps = rec.tid == it->tid && rec.pos0 < it->end && rec.end0 > it->beg;
        if (overlaps && !seen_has(it, rec.aln_id)) {
            seen_add(it, rec.aln_id);
            *out = rec;
            return 1;
        }
        lmdb_view_record_clear(&rec);
        if ((it->bin + 1) * bin_size < it->beg) {
            ++it->bin;
            it->started = 0;
        }
    }
    return 0;
}

int lmdb_view_iter_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    memset(out, 0, sizeof *out);
    if (it->kind == LMDB_VIEW_ITER_NAME) return iter_name_next(it, out);
    if (it->kind == LMDB_VIEW_ITER_COORD) return iter_coord_next(it, out);
    if (it->kind == LMDB_VIEW_ITER_REGION) return iter_region_next(it, out);
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    if (!iter_cursor_next(it, &key, &val)) return 0;
    return iter_decode_ref(it, &val, out, "aln_import_by_order");
}

void lmdb_view_iter_close(lmdb_view_iter_t *it) {
    if (!it) return;
    if (it->cursor) mdb_cursor_close(it->cursor);
    if (it->read_cursor) mdb_cursor_close(it->read_cursor);
    if (it->txn) mdb_txn_abort(it->txn);
    free(it->seen_aln_ids);
    free(it);
}

void lmdb_view_record_clear(lmdb_view_record_t *rec) {
    if (!rec) return;
    free(rec->cigar);
    free(rec->aux);
    memset(rec, 0, sizeof *rec);
}

#endif
