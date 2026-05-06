#include "common.h"
#include "backend_lmdb.h"

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

static MDB_val subset_read_key(int64_t subset_id, int64_t read_id) {
    unsigned char *buf = xmalloc(16);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)subset_id);
    put_u64(&p, (uint64_t)read_id);
    MDB_val key = { 16, buf };
    return key;
}

static void subset_read_key_parts(const MDB_val *key,
                                  int64_t *subset_id,
                                  int64_t *read_id) {
    if (key->mv_size < 16) die("corrupt lmdb subset_read key");
    const unsigned char *p = key->mv_data;
    *subset_id = (int64_t)get_u64(&p);
    *read_id = (int64_t)get_u64(&p);
}

static MDB_val read_tag_key(int64_t read_id, const char *tag) {
    if (strlen(tag) != 2) die("invalid read tag '%s'", tag);
    unsigned char *buf = xmalloc(10);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)read_id);
    *p++ = (unsigned char)tag[0];
    *p++ = (unsigned char)tag[1];
    MDB_val key = { 10, buf };
    return key;
}

static void read_tag_key_parts(const MDB_val *key,
                               int64_t *read_id,
                               char tag[3]) {
    if (key->mv_size < 10) die("corrupt lmdb read_tag key");
    const unsigned char *p = key->mv_data;
    *read_id = (int64_t)get_u64(&p);
    tag[0] = (char)*p++;
    tag[1] = (char)*p++;
    tag[2] = 0;
}

static MDB_val checkpoint_key(int tid, int64_t checkpoint_pos) {
    unsigned char *buf = xmalloc(16);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)tid);
    put_u64(&p, (uint64_t)checkpoint_pos);
    MDB_val key = { 16, buf };
    return key;
}

static void checkpoint_key_parts(const MDB_val *key,
                                 int *tid,
                                 int64_t *checkpoint_pos) {
    if (key->mv_size < 16) die("corrupt lmdb checkpoint key");
    const unsigned char *p = key->mv_data;
    *tid = (int32_t)get_u64(&p);
    *checkpoint_pos = (int64_t)get_u64(&p);
}

static MDB_val i64_value(int64_t value) {
    unsigned char *buf = xmalloc(8);
    unsigned char *p = buf;
    put_u64(&p, (uint64_t)value);
    MDB_val val = { 8, buf };
    return val;
}

static MDB_val empty_value(void) {
    MDB_val val = { 0, NULL };
    return val;
}

static void coord_key_parts(const MDB_val *key, int *tid, int64_t *pos0) {
    if (key->mv_size < 16) die("corrupt lmdb coord key");
    const unsigned char *p = key->mv_data;
    *tid = (int32_t)get_u64(&p);
    *pos0 = (int64_t)get_u64(&p);
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

static void decode_view_read(const MDB_val *val,
                             int64_t *store_read_id,
                             char **qname,
                             int64_t *hash) {
    const unsigned char *p = val->mv_data;
    const unsigned char *end = p + val->mv_size;
    if (p >= end || *p++ != 1) die("unsupported lmdb view read version");
    uint16_t qn = get_u16(&p);
    if (p + qn + 16 > end) die("truncated lmdb view read record");
    *qname = xmalloc((size_t)qn + 1);
    memcpy(*qname, p, qn);
    (*qname)[qn] = 0;
    p += qn;
    *store_read_id = (int64_t)get_u64(&p);
    *hash = (int64_t)get_u64(&p);
}

static MDB_val encode_subset_value(int64_t subset_id,
                                   const char *description,
                                   const char *created_at) {
    size_t dn = description ? strlen(description) : 0;
    size_t cn = created_at ? strlen(created_at) : 0;
    if (dn > UINT16_MAX || cn > UINT16_MAX)
        die("subset metadata value too long");
    size_t n = 1 + 8 + 2 + dn + 2 + cn;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    put_u64(&p, (uint64_t)subset_id);
    put_u16(&p, (uint16_t)dn);
    if (dn) { memcpy(p, description, dn); p += dn; }
    put_u16(&p, (uint16_t)cn);
    if (cn) memcpy(p, created_at, cn);
    MDB_val val = { n, buf };
    return val;
}

static int64_t decode_subset_id(const MDB_val *val) {
    const unsigned char *p = val->mv_data;
    const unsigned char *end = p + val->mv_size;
    if (p >= end || *p++ != 1 || p + 8 > end)
        die("corrupt lmdb subset record");
    return (int64_t)get_u64(&p);
}

static MDB_val encode_read_tag_value(const char *type,
                                     const void *value,
                                     size_t value_len,
                                     const char *source) {
    if (!type || strlen(type) != 1) die("invalid read tag type");
    size_t sn = source ? strlen(source) : 0;
    if (value_len > UINT32_MAX || sn > UINT16_MAX)
        die("read tag value too long");
    size_t n = 1 + 1 + 4 + value_len + 2 + sn;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    *p++ = (unsigned char)type[0];
    put_u32(&p, (uint32_t)value_len);
    if (value_len) { memcpy(p, value, value_len); p += value_len; }
    put_u16(&p, (uint16_t)sn);
    if (sn) memcpy(p, source, sn);
    MDB_val val = { n, buf };
    return val;
}

static MDB_val encode_annotation_meta_value(const char *kind,
                                            const char *source,
                                            const char *command_line,
                                            const char *created_at) {
    size_t kn = kind ? strlen(kind) : 0;
    size_t sn = source ? strlen(source) : 0;
    size_t cn = command_line ? strlen(command_line) : 0;
    size_t tn = created_at ? strlen(created_at) : 0;
    if (kn > UINT16_MAX || sn > UINT16_MAX || cn > UINT16_MAX || tn > UINT16_MAX)
        die("annotation metadata value too long");
    size_t n = 1 + 2 + kn + 2 + sn + 2 + cn + 2 + tn;
    unsigned char *buf = xmalloc(n);
    unsigned char *p = buf;
    *p++ = 1;
    put_u16(&p, (uint16_t)kn); if (kn) { memcpy(p, kind, kn); p += kn; }
    put_u16(&p, (uint16_t)sn); if (sn) { memcpy(p, source, sn); p += sn; }
    put_u16(&p, (uint16_t)cn); if (cn) { memcpy(p, command_line, cn); p += cn; }
    put_u16(&p, (uint16_t)tn); if (tn) memcpy(p, created_at, tn);
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
    dbi_open_optional(txn, path, "aln_checkpoint", dbi_flags,
                      &b->aln_checkpoint);
    dbi_open_optional(txn, path, "subset_by_name", dbi_flags,
                      &b->subset_by_name);
    dbi_open_optional(txn, path, "subset_read", dbi_flags,
                      &b->subset_read);
    dbi_open_optional(txn, path, "read_tag", dbi_flags, &b->read_tag);
    dbi_open_optional(txn, path, "annotation_meta", dbi_flags,
                      &b->annotation_meta);

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
        if (b->aln_checkpoint) mdb_dbi_close(b->env, b->aln_checkpoint);
        if (b->subset_by_name) mdb_dbi_close(b->env, b->subset_by_name);
        if (b->subset_read) mdb_dbi_close(b->env, b->subset_read);
        if (b->read_tag) mdb_dbi_close(b->env, b->read_tag);
        if (b->annotation_meta) mdb_dbi_close(b->env, b->annotation_meta);
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

static char *lmdb_backend_meta_get_txn(lmdb_backend_t *b,
                                       MDB_txn *txn,
                                       const char *key) {
    MDB_val k = { strlen(key), (void *)key };
    MDB_val v = { 0, NULL };
    int rc = mdb_get(txn, b->meta, &k, &v);
    if (rc == MDB_NOTFOUND) return NULL;
    if (rc) lmdb_die(rc, "get(meta)", b->path);

    char *out = xmalloc(v.mv_size + 1);
    memcpy(out, v.mv_data, v.mv_size);
    out[v.mv_size] = 0;
    return out;
}

static int lmdb_backend_meta_get_i64_txn(lmdb_backend_t *b,
                                         MDB_txn *txn,
                                         const char *key,
                                         int64_t *out) {
    char *s = lmdb_backend_meta_get_txn(b, txn, key);
    if (!s) return 0;
    *out = (int64_t)strtoll(s, NULL, 10);
    free(s);
    return 1;
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
    if (clear_interval) {
        lmdb_drop_clear(b, txn, b->aln_by_bin, "drop(aln_by_bin)");
        lmdb_drop_clear(b, txn, b->aln_checkpoint, "drop(aln_checkpoint)");
    }
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

    (void)index_interval;

    write_txn_finish(b, txn, owned);
}

typedef struct active_interval {
    int64_t start;
    int64_t end;
} active_interval_t;

static int64_t active_min_start(active_interval_t *active,
                                size_t *n_active,
                                int64_t checkpoint_pos) {
    int64_t min_start = checkpoint_pos;
    size_t out = 0;
    int have_active = 0;
    for (size_t i = 0; i < *n_active; ++i) {
        if (active[i].end <= checkpoint_pos) continue;
        if (!have_active || active[i].start < min_start) {
            min_start = active[i].start;
            have_active = 1;
        }
        active[out++] = active[i];
    }
    *n_active = out;
    return min_start;
}

static void checkpoint_put(lmdb_backend_t *b, MDB_txn *txn,
                           int tid, int64_t checkpoint_pos,
                           int64_t min_active_start) {
    MDB_val key = checkpoint_key(tid, checkpoint_pos);
    MDB_val val = i64_value(min_active_start);
    int rc = mdb_put(txn, b->aln_checkpoint, &key, &val, 0);
    if (rc) lmdb_die(rc, "put(aln_checkpoint)", b->path);
    free(key.mv_data);
    free(val.mv_data);
}

void lmdb_view_build_checkpoints(lmdb_backend_t *b, int64_t checkpoint_step) {
    if (checkpoint_step <= 0) die("invalid LMDB checkpoint step");

    MDB_txn *read_txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &read_txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    MDB_cursor *cursor = NULL;
    rc = mdb_cursor_open(read_txn, b->aln_by_coord, &cursor);
    if (rc) lmdb_die(rc, "cursor(aln_by_coord)", b->path);

    lmdb_backend_begin(b);
    MDB_txn *write = b->write_txn;

    active_interval_t *active = NULL;
    size_t n_active = 0, cap_active = 0;
    int current_tid = -1;
    int64_t next_checkpoint = 0;
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while ((rc = mdb_cursor_get(cursor, &key, &val,
                               key.mv_data ? MDB_NEXT : MDB_FIRST)) == 0) {
        int tid;
        int64_t pos0;
        coord_key_parts(&key, &tid, &pos0);
        int64_t aln_id = decode_i64_value(&val, "aln_by_coord");
        lmdb_view_record_t rec;
        if (!lmdb_view_get_record_txn(b, read_txn, aln_id, &rec))
            die("lmdb alignment record '%lld' referenced by aln_by_coord is missing in '%s'",
                (long long)aln_id, b->path);

        if (tid != current_tid) {
            current_tid = tid;
            n_active = 0;
            next_checkpoint = (pos0 / checkpoint_step) * checkpoint_step;
        }
        while (next_checkpoint <= pos0) {
            int64_t min_start = active_min_start(active, &n_active,
                                                 next_checkpoint);
            checkpoint_put(b, write, current_tid, next_checkpoint, min_start);
            next_checkpoint += checkpoint_step;
        }
        if (rec.end0 > rec.pos0) {
            if (n_active == cap_active) {
                cap_active = cap_active ? cap_active * 2 : 64;
                active = xrealloc(active, cap_active * sizeof *active);
            }
            active[n_active++] = (active_interval_t){ rec.pos0, rec.end0 };
        }
        lmdb_view_record_clear(&rec);
    }
    if (rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(aln_by_coord)", b->path);

    free(active);
    lmdb_backend_commit(b);
    mdb_cursor_close(cursor);
    mdb_txn_abort(read_txn);
}

int64_t lmdb_subset_upsert(lmdb_backend_t *b,
                           const char *name,
                           const char *description) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    MDB_val key = { strlen(name), (void *)name };
    MDB_val old = { 0, NULL };
    int rc = mdb_get(txn, b->subset_by_name, &key, &old);
    int64_t subset_id;
    if (rc == MDB_NOTFOUND) {
        subset_id = 1;
        char *next = lmdb_backend_meta_get_txn(b, txn, "subset_next_id");
        if (next) {
            subset_id = (int64_t)strtoll(next, NULL, 10);
            free(next);
        }
        char buf[32];
        snprintf(buf, sizeof buf, "%lld", (long long)(subset_id + 1));
        MDB_val nk = { strlen("subset_next_id"), "subset_next_id" };
        MDB_val nv = { strlen(buf), buf };
        rc = mdb_put(txn, b->meta, &nk, &nv, 0);
        if (rc) lmdb_die(rc, "put(meta)", b->path);
    } else if (rc == 0) {
        subset_id = decode_subset_id(&old);
    } else {
        lmdb_die(rc, "get(subset_by_name)", b->path);
    }

    char *ts = now_iso8601();
    MDB_val val = encode_subset_value(subset_id, description, ts);
    free(ts);
    rc = mdb_put(txn, b->subset_by_name, &key, &val, 0);
    if (rc) lmdb_die(rc, "put(subset_by_name)", b->path);
    free(val.mv_data);
    write_txn_finish(b, txn, owned);
    return subset_id;
}

int lmdb_subset_lookup(lmdb_backend_t *b,
                       const char *name,
                       int64_t *subset_id) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    MDB_val key = { strlen(name), (void *)name };
    MDB_val val = { 0, NULL };
    rc = mdb_get(txn, b->subset_by_name, &key, &val);
    if (rc == MDB_NOTFOUND) {
        mdb_txn_abort(txn);
        return 0;
    }
    if (rc) lmdb_die(rc, "get(subset_by_name)", b->path);
    *subset_id = decode_subset_id(&val);
    mdb_txn_abort(txn);
    return 1;
}

void lmdb_subset_add_read(lmdb_backend_t *b,
                          int64_t subset_id,
                          int64_t read_id) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    MDB_val key = subset_read_key(subset_id, read_id);
    MDB_val val = empty_value();
    int rc = mdb_put(txn, b->subset_read, &key, &val, MDB_NOOVERWRITE);
    if (rc && rc != MDB_KEYEXIST) lmdb_die(rc, "put(subset_read)", b->path);
    free(key.mv_data);
    write_txn_finish(b, txn, owned);
}

int lmdb_subset_contains_read(lmdb_backend_t *b,
                              int64_t subset_id,
                              int64_t read_id) {
    MDB_txn *txn = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    MDB_val key = subset_read_key(subset_id, read_id);
    MDB_val val = { 0, NULL };
    rc = mdb_get(txn, b->subset_read, &key, &val);
    free(key.mv_data);
    mdb_txn_abort(txn);
    if (rc == MDB_NOTFOUND) return 0;
    if (rc) lmdb_die(rc, "get(subset_read)", b->path);
    return 1;
}

void lmdb_read_tag_put(lmdb_backend_t *b,
                       int64_t read_id,
                       const char *tag,
                       const char *type,
                       const void *value,
                       size_t value_len,
                       const char *source) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    MDB_val key = read_tag_key(read_id, tag);
    MDB_val val = encode_read_tag_value(type, value, value_len, source);
    int rc = mdb_put(txn, b->read_tag, &key, &val, 0);
    if (rc) lmdb_die(rc, "put(read_tag)", b->path);
    free(key.mv_data);
    free(val.mv_data);
    write_txn_finish(b, txn, owned);
}

void lmdb_annotation_meta_put(lmdb_backend_t *b,
                              const char *kind,
                              const char *source,
                              const char *command_line) {
    int owned = 0;
    MDB_txn *txn = write_txn(b, &owned);
    int64_t annotation_id = 1;
    char *next = lmdb_backend_meta_get_txn(b, txn, "annotation_next_id");
    if (next) {
        annotation_id = (int64_t)strtoll(next, NULL, 10);
        free(next);
    }
    char nextbuf[32];
    snprintf(nextbuf, sizeof nextbuf, "%lld", (long long)(annotation_id + 1));
    MDB_val mk = { strlen("annotation_next_id"), "annotation_next_id" };
    MDB_val mv = { strlen(nextbuf), nextbuf };
    int rc = mdb_put(txn, b->meta, &mk, &mv, 0);
    if (rc) lmdb_die(rc, "put(meta)", b->path);

    unsigned char keybuf[8];
    be64_key(annotation_id, keybuf);
    MDB_val key = { sizeof keybuf, keybuf };
    char *ts = now_iso8601();
    MDB_val val = encode_annotation_meta_value(kind, source, command_line, ts);
    free(ts);
    rc = mdb_put(txn, b->annotation_meta, &key, &val, 0);
    if (rc) lmdb_die(rc, "put(annotation_meta)", b->path);
    free(val.mv_data);
    write_txn_finish(b, txn, owned);
}

typedef struct rehint_map_entry {
    int64_t old_id;
    int64_t new_id;
    int64_t qname_hash;
    char *qname;
} rehint_map_entry_t;

typedef struct subset_move_entry {
    int64_t subset_id;
    int64_t read_id;
} subset_move_entry_t;

typedef struct tag_move_entry {
    int64_t read_id;
    char tag[3];
    void *value;
    size_t value_len;
} tag_move_entry_t;

static int64_t rehint_map_lookup(const rehint_map_entry_t *map,
                                 size_t n_map,
                                 int64_t old_id) {
    size_t lo = 0, hi = n_map;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (map[mid].old_id < old_id) lo = mid + 1;
        else hi = mid;
    }
    if (lo == n_map || map[lo].old_id != old_id)
        die("rehint: store_read_id '%lld' is absent from read_dict",
            (long long)old_id);
    return map[lo].new_id;
}

static int64_t lookup_store_read_id_txn(lmdb_backend_t *store,
                                        MDB_txn *txn,
                                        int64_t hash,
                                        const char *qname) {
    MDB_val key = qname_key(hash, qname);
    MDB_val val = { 0, NULL };
    int rc = mdb_get(txn, store->canonical_by_qname, &key, &val);
    free(key.mv_data);
    if (rc == MDB_NOTFOUND)
        die("rehint: qname '%s' is absent from store index '%s'",
            qname, store->path);
    if (rc) lmdb_die(rc, "get(canonical_by_qname)", store->path);
    const unsigned char *p = val.mv_data;
    return (int64_t)get_u64(&p);
}

void lmdb_view_rehint(lmdb_backend_t *view,
                      lmdb_backend_t *store,
                      const char *store_path,
                      const char *store_index_path,
                      int64_t *n_reads,
                      int64_t *n_changed) {
    if (n_reads) *n_reads = 0;
    if (n_changed) *n_changed = 0;

    MDB_txn *store_txn = NULL;
    int rc = mdb_txn_begin(store->env, NULL, MDB_RDONLY, &store_txn);
    if (rc) lmdb_die(rc, "txn_begin", store->path);

    lmdb_backend_begin(view);
    MDB_txn *txn = view->write_txn;

    rehint_map_entry_t *map = NULL;
    size_t n_map = 0, cap_map = 0;
    MDB_cursor *cursor = NULL;
    rc = mdb_cursor_open(txn, view->read_by_id, &cursor);
    if (rc) lmdb_die(rc, "cursor(read_by_id)", view->path);
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while ((rc = mdb_cursor_get(cursor, &key, &val,
                               n_map ? MDB_NEXT : MDB_FIRST)) == 0) {
        int64_t old_id, stored_id, hash;
        char *qname = NULL;
        const unsigned char *kp = key.mv_data;
        old_id = (int64_t)get_u64(&kp);
        decode_view_read(&val, &stored_id, &qname, &hash);
        if (stored_id != old_id)
            die("rehint: corrupt read_by_id entry in '%s'", view->path);
        int64_t new_id = lookup_store_read_id_txn(store, store_txn, hash, qname);
        if (n_map == cap_map) {
            cap_map = cap_map ? cap_map * 2 : 1024;
            map = xrealloc(map, cap_map * sizeof *map);
        }
        map[n_map++] = (rehint_map_entry_t){ old_id, new_id, hash, qname };
        if (old_id != new_id && n_changed) ++*n_changed;
    }
    if (rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(read_by_id)", view->path);
    mdb_cursor_close(cursor);
    if (n_reads) *n_reads = (int64_t)n_map;

    subset_move_entry_t *subset_moves = NULL;
    size_t n_subset = 0, cap_subset = 0;
    rc = mdb_cursor_open(txn, view->subset_read, &cursor);
    if (rc) lmdb_die(rc, "cursor(subset_read)", view->path);
    int started = 0;
    while ((rc = mdb_cursor_get(cursor, &key, &val,
                               started ? MDB_NEXT : MDB_FIRST)) == 0) {
        started = 1;
        int64_t subset_id, old_id;
        subset_read_key_parts(&key, &subset_id, &old_id);
        if (n_subset == cap_subset) {
            cap_subset = cap_subset ? cap_subset * 2 : 1024;
            subset_moves = xrealloc(subset_moves, cap_subset * sizeof *subset_moves);
        }
        subset_moves[n_subset++] = (subset_move_entry_t){
            subset_id, rehint_map_lookup(map, n_map, old_id)
        };
    }
    if (rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(subset_read)", view->path);
    mdb_cursor_close(cursor);

    tag_move_entry_t *tag_moves = NULL;
    size_t n_tag = 0, cap_tag = 0;
    rc = mdb_cursor_open(txn, view->read_tag, &cursor);
    if (rc) lmdb_die(rc, "cursor(read_tag)", view->path);
    started = 0;
    while ((rc = mdb_cursor_get(cursor, &key, &val,
                               started ? MDB_NEXT : MDB_FIRST)) == 0) {
        started = 1;
        int64_t old_id;
        char tag[3];
        read_tag_key_parts(&key, &old_id, tag);
        if (n_tag == cap_tag) {
            cap_tag = cap_tag ? cap_tag * 2 : 1024;
            tag_moves = xrealloc(tag_moves, cap_tag * sizeof *tag_moves);
        }
        tag_moves[n_tag].read_id = rehint_map_lookup(map, n_map, old_id);
        memcpy(tag_moves[n_tag].tag, tag, sizeof tag_moves[n_tag].tag);
        tag_moves[n_tag].value_len = val.mv_size;
        tag_moves[n_tag].value = xmalloc(val.mv_size);
        memcpy(tag_moves[n_tag].value, val.mv_data, val.mv_size);
        ++n_tag;
    }
    if (rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(read_tag)", view->path);
    mdb_cursor_close(cursor);

    lmdb_drop_clear(view, txn, view->read_by_id, "drop(read_by_id)");
    lmdb_drop_clear(view, txn, view->read_by_qname, "drop(read_by_qname)");
    lmdb_drop_clear(view, txn, view->aln_by_read, "drop(aln_by_read)");
    lmdb_drop_clear(view, txn, view->subset_read, "drop(subset_read)");
    lmdb_drop_clear(view, txn, view->read_tag, "drop(read_tag)");

    for (size_t i = 0; i < n_map; ++i)
        lmdb_view_put_read(view, map[i].new_id, map[i].qname, map[i].qname_hash);

    rc = mdb_cursor_open(txn, view->aln_by_id, &cursor);
    if (rc) lmdb_die(rc, "cursor(aln_by_id)", view->path);
    started = 0;
    while ((rc = mdb_cursor_get(cursor, &key, &val,
                               started ? MDB_NEXT : MDB_FIRST)) == 0) {
        started = 1;
        lmdb_view_record_t rec;
        decode_view_record(&val, &rec);
        rec.store_read_id = rehint_map_lookup(map, n_map, rec.store_read_id);
        MDB_val new_val = encode_view_record(&rec);
        rc = mdb_cursor_put(cursor, &key, &new_val, MDB_CURRENT);
        if (rc) lmdb_die(rc, "put(aln_by_id)", view->path);

        unsigned char refbuf[8];
        unsigned char *p = refbuf;
        put_u64(&p, (uint64_t)rec.aln_id);
        MDB_val refval = { sizeof refbuf, refbuf };
        MDB_val rkey = read_aln_key(rec.store_read_id, rec.record_order,
                                    rec.aln_id);
        rc = mdb_put(txn, view->aln_by_read, &rkey, &refval, 0);
        if (rc) lmdb_die(rc, "put(aln_by_read)", view->path);
        free(rkey.mv_data);
        free(new_val.mv_data);
        lmdb_view_record_clear(&rec);
    }
    if (rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(aln_by_id)", view->path);
    mdb_cursor_close(cursor);

    for (size_t i = 0; i < n_subset; ++i) {
        MDB_val skey = subset_read_key(subset_moves[i].subset_id,
                                       subset_moves[i].read_id);
        MDB_val empty = empty_value();
        rc = mdb_put(txn, view->subset_read, &skey, &empty, MDB_NOOVERWRITE);
        if (rc && rc != MDB_KEYEXIST) lmdb_die(rc, "put(subset_read)", view->path);
        free(skey.mv_data);
    }

    for (size_t i = 0; i < n_tag; ++i) {
        MDB_val tkey = read_tag_key(tag_moves[i].read_id, tag_moves[i].tag);
        MDB_val tval = { tag_moves[i].value_len, tag_moves[i].value };
        rc = mdb_put(txn, view->read_tag, &tkey, &tval, 0);
        if (rc) lmdb_die(rc, "put(read_tag)", view->path);
        free(tkey.mv_data);
    }

    if (store_path) lmdb_backend_meta_set(view, "store_path", store_path);
    if (store_index_path) lmdb_backend_meta_set(view, "store_index_path", store_index_path);
    char *shh = lmdb_backend_meta_get_txn(store, store_txn, "store_header_hash");
    if (shh) { lmdb_backend_meta_set(view, "store_header_hash", shh); free(shh); }
    char *ssh = lmdb_backend_meta_get_txn(store, store_txn, "store_sq_hash");
    if (ssh) { lmdb_backend_meta_set(view, "store_sq_hash", ssh); free(ssh); }
    char *ts = now_iso8601();
    lmdb_backend_meta_set(view, "rehinted_at", ts);
    free(ts);
    lmdb_view_put_meta_i64(view, "rehinted_read_count", (int64_t)n_map);

    lmdb_backend_commit(view);
    mdb_txn_abort(store_txn);

    for (size_t i = 0; i < n_map; ++i) free(map[i].qname);
    for (size_t i = 0; i < n_tag; ++i) free(tag_moves[i].value);
    free(map);
    free(subset_moves);
    free(tag_moves);
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
    int64_t scan_start;
    unsigned char read_prefix[8];
    int have_read_prefix;
    unsigned char subset_prefix[8];
    int have_subset_prefix;
};

static int64_t max_span_scan_start(lmdb_backend_t *b,
                                   MDB_txn *txn,
                                   int64_t beg) {
    int64_t max_span = 0;
    lmdb_backend_meta_get_i64_txn(b, txn, "view_max_ref_span", &max_span);
    if (max_span < 0) max_span = 0;
    return beg > max_span ? beg - max_span : 0;
}

static int64_t checkpoint_scan_start(lmdb_backend_t *b,
                                     MDB_txn *txn,
                                     int tid,
                                     int64_t beg) {
    char *strategy = lmdb_backend_meta_get_txn(b, txn, "interval_strategy");
    int use_checkpoint = strategy && strcmp(strategy, "checkpoint") == 0;
    free(strategy);
    if (!use_checkpoint) return max_span_scan_start(b, txn, beg);

    MDB_cursor *cursor = NULL;
    int rc = mdb_cursor_open(txn, b->aln_checkpoint, &cursor);
    if (rc) lmdb_die(rc, "cursor(aln_checkpoint)", b->path);

    unsigned char searchbuf[16];
    unsigned char *sp = searchbuf;
    put_u64(&sp, (uint64_t)tid);
    put_u64(&sp, (uint64_t)(beg + 1));
    MDB_val key = { sizeof searchbuf, searchbuf };
    MDB_val val = { 0, NULL };
    rc = mdb_cursor_get(cursor, &key, &val, MDB_SET_RANGE);
    if (rc == MDB_NOTFOUND)
        rc = mdb_cursor_get(cursor, &key, &val, MDB_LAST);
    else if (rc == 0)
        rc = mdb_cursor_get(cursor, &key, &val, MDB_PREV);
    if (rc == MDB_NOTFOUND) {
        mdb_cursor_close(cursor);
        return 0;
    }
    if (rc) lmdb_die(rc, "cursor(aln_checkpoint)", b->path);

    int key_tid;
    int64_t checkpoint_pos;
    checkpoint_key_parts(&key, &key_tid, &checkpoint_pos);
    if (key_tid != tid || checkpoint_pos > beg) {
        mdb_cursor_close(cursor);
        return max_span_scan_start(b, txn, beg);
    }
    int64_t scan_start = decode_i64_value(&val, "aln_checkpoint");
    mdb_cursor_close(cursor);
    return scan_start < 0 ? 0 : scan_start;
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
        dbi = b->aln_by_coord;
        it->scan_start = checkpoint_scan_start(b, it->txn, tid, beg);
    } else if (kind == LMDB_VIEW_ITER_READ_ID) {
        dbi = b->aln_by_read;
    } else if (kind == LMDB_VIEW_ITER_NAME) {
        rc = mdb_cursor_open(it->txn, b->read_by_qname, &it->read_cursor);
        if (rc) lmdb_die(rc, "cursor(read_by_qname)", b->path);
        dbi = b->aln_by_read;
    }
    rc = mdb_cursor_open(it->txn, dbi, &it->cursor);
    if (rc) lmdb_die(rc, "cursor(view)", b->path);
    return it;
}

lmdb_view_iter_t *lmdb_view_iter_open_read_id(lmdb_backend_t *b,
                                              int64_t store_read_id) {
    lmdb_view_iter_t *it = lmdb_view_iter_open(b, LMDB_VIEW_ITER_READ_ID,
                                               -1, 0, 0);
    unsigned char *p = it->read_prefix;
    put_u64(&p, (uint64_t)store_read_id);
    it->have_read_prefix = 1;
    return it;
}

lmdb_view_iter_t *lmdb_view_iter_open_subset(lmdb_backend_t *b,
                                             int64_t subset_id) {
    lmdb_view_iter_t *it = xcalloc(1, sizeof *it);
    it->b = b;
    it->kind = LMDB_VIEW_ITER_SUBSET;
    unsigned char *p = it->subset_prefix;
    put_u64(&p, (uint64_t)subset_id);
    it->have_subset_prefix = 1;

    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &it->txn);
    if (rc) lmdb_die(rc, "txn_begin", b->path);
    rc = mdb_cursor_open(it->txn, b->subset_read, &it->read_cursor);
    if (rc) lmdb_die(rc, "cursor(subset_read)", b->path);
    rc = mdb_cursor_open(it->txn, b->aln_by_read, &it->cursor);
    if (rc) lmdb_die(rc, "cursor(aln_by_read)", b->path);
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
            if (!it->started) {
                key.mv_size = sizeof it->read_prefix;
                key.mv_data = it->read_prefix;
            }
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

static int iter_read_id_next(lmdb_view_iter_t *it,
                             lmdb_view_record_t *out) {
    if (!it->have_read_prefix) return 0;
    MDB_val key = { sizeof it->read_prefix, it->read_prefix };
    MDB_val val = { 0, NULL };
    int rc = mdb_cursor_get(it->cursor, &key, &val,
                            it->started ? MDB_NEXT : MDB_SET_RANGE);
    it->started = 1;
    if (rc == MDB_NOTFOUND) return 0;
    if (rc) lmdb_die(rc, "cursor(aln_by_read)", it->b->path);
    if (!key_has_prefix(&key, it->read_prefix, sizeof it->read_prefix))
        return 0;
    return iter_decode_ref(it, &val, out, "aln_by_read");
}

static int iter_subset_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while (1) {
        if (it->have_read_prefix) {
            if (!it->started) {
                key.mv_size = sizeof it->read_prefix;
                key.mv_data = it->read_prefix;
            }
            int rc = mdb_cursor_get(it->cursor, &key, &val,
                                    it->started ? MDB_NEXT : MDB_SET_RANGE);
            it->started = 1;
            if (rc && rc != MDB_NOTFOUND) lmdb_die(rc, "cursor(aln_by_read)", it->b->path);
            if (!rc && key_has_prefix(&key, it->read_prefix, sizeof it->read_prefix))
                return iter_decode_ref(it, &val, out, "aln_by_read");
            it->have_read_prefix = 0;
            it->started = 0;
        }

        int rc;
        if (it->phase) {
            rc = mdb_cursor_get(it->read_cursor, &key, &val, MDB_NEXT);
        } else {
            key.mv_size = sizeof it->subset_prefix;
            key.mv_data = it->subset_prefix;
            rc = mdb_cursor_get(it->read_cursor, &key, &val, MDB_SET_RANGE);
        }
        it->phase = 1;
        if (rc == MDB_NOTFOUND) return 0;
        if (rc) lmdb_die(rc, "cursor(subset_read)", it->b->path);
        if (!key_has_prefix(&key, it->subset_prefix, sizeof it->subset_prefix))
            return 0;

        int64_t subset_id, read_id;
        subset_read_key_parts(&key, &subset_id, &read_id);
        (void)subset_id;
        unsigned char *rp = it->read_prefix;
        put_u64(&rp, (uint64_t)read_id);
        it->have_read_prefix = 1;
        it->started = 0;
    }
}

static int iter_region_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    MDB_val key = { 0, NULL }, val = { 0, NULL };
    while (1) {
        int rc;
        if (!it->started) {
            unsigned char searchbuf[32];
            unsigned char *sp = searchbuf;
            put_u64(&sp, (uint64_t)it->tid);
            put_u64(&sp, (uint64_t)it->scan_start);
            put_u64(&sp, 0);
            put_u64(&sp, 0);
            key.mv_size = sizeof searchbuf;
            key.mv_data = searchbuf;
            rc = mdb_cursor_get(it->cursor, &key, &val, MDB_SET_RANGE);
            it->started = 1;
        } else {
            rc = mdb_cursor_get(it->cursor, &key, &val, MDB_NEXT);
        }
        if (rc == MDB_NOTFOUND) return 0;
        if (rc) lmdb_die(rc, "cursor(aln_by_coord)", it->b->path);

        int key_tid;
        int64_t pos0;
        coord_key_parts(&key, &key_tid, &pos0);
        if (key_tid != it->tid || pos0 >= it->end) return 0;

        int64_t aln_id = decode_i64_value(&val, "aln_by_coord");
        lmdb_view_record_t rec;
        if (!lmdb_view_get_record_txn(it->b, it->txn, aln_id, &rec))
            die("lmdb alignment record '%lld' referenced by aln_by_coord is missing in '%s'",
                (long long)aln_id, it->b->path);
        int overlaps = rec.tid == it->tid && rec.pos0 < it->end && rec.end0 > it->beg;
        if (overlaps) {
            *out = rec;
            return 1;
        }
        lmdb_view_record_clear(&rec);
    }
}

int lmdb_view_iter_next(lmdb_view_iter_t *it, lmdb_view_record_t *out) {
    memset(out, 0, sizeof *out);
    if (it->kind == LMDB_VIEW_ITER_NAME) return iter_name_next(it, out);
    if (it->kind == LMDB_VIEW_ITER_READ_ID) return iter_read_id_next(it, out);
    if (it->kind == LMDB_VIEW_ITER_SUBSET) return iter_subset_next(it, out);
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
    free(it);
}

void lmdb_view_record_clear(lmdb_view_record_t *rec) {
    if (!rec) return;
    free(rec->cigar);
    free(rec->aux);
    memset(rec, 0, sizeof *rec);
}
