#include "store.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <errno.h>

typedef struct payload_cache_entry {
    int64_t  read_id;
    int      next;
    bam1_t  *bam;
} payload_cache_entry_t;

struct payload_cache {
    payload_cache_entry_t *entries;
    int *buckets;
    size_t n_entries;
    size_t cap_entries;
    size_t bucket_mask;
    size_t evict_cursor;
    size_t max_bytes;
    size_t used_bytes;
};

static size_t bam_payload_bytes(const bam1_t *b) {
    return b ? sizeof *b + (size_t)b->m_data : 0;
}

store_t *store_open(const char *store_path, const char *index_path) {
    store_t *st = xcalloc(1, sizeof *st);
    st->store_path = xstrdup(store_path);
    st->index_path = xstrdup(index_path);

    struct stat idx_st;
    if (stat(index_path, &idx_st) != 0)
        die("stat store index '%s' failed: %s", index_path, strerror(errno));
    if (S_ISDIR(idx_st.st_mode))
        die("store index '%s' is a directory; LMDB store indexes are file-backed .ssam files",
            index_path);
    st->lmdb = lmdb_backend_open(index_path, /*readonly*/1, 0);
    char *type = lmdb_backend_meta_get(st->lmdb, "isam_type");
    if (!type || strcmp(type, "store") != 0)
        die("'%s' is not an LMDB store index (type=%s)",
            index_path, type ? type : "<missing>");
    free(type);

    st->sam = sam_open(store_path, "r");
    if (!st->sam) die("cannot open canonical store '%s'", store_path);
    const htsFormat *fmt = hts_get_format(st->sam);
    st->is_bam = (fmt && fmt->format == bam);

    st->hdr = sam_hdr_read(st->sam);
    if (!st->hdr) die("cannot read header of canonical store '%s'", store_path);

    return st;
}

void store_close(store_t *st) {
    if (!st) return;
    if (st->lmdb) lmdb_backend_close(st->lmdb);
    if (st->hdr)      sam_hdr_destroy(st->hdr);
    if (st->sam)      sam_close(st->sam);
    free(st->store_path);
    free(st->index_path);
    free(st);
}

int store_has_name_index(store_t *st) {
    int64_t indexed_name = 0;
    return lmdb_backend_meta_get_i64(st->lmdb, "indexed_name", &indexed_name) &&
           indexed_name;
}

static void fill_record_from_lmdb(lmdb_store_record_t *lr, store_record_t *out) {
    memset(out, 0, sizeof *out);
    out->read_id = lr->read_id;
    out->qname = xstrdup(lr->qname);
    out->tid = lr->tid;
    out->pos0 = lr->pos0;
    out->end0 = lr->end0;
    out->flag = lr->flag;
    out->is_unmapped = lr->is_unmapped;
    out->bam_voffset = lr->bam_voffset;
    out->has_voffset = lr->has_voffset;
}

int store_lookup_qname(store_t *st, const char *qname, store_record_t *out) {
    lmdb_store_record_t lr;
    if (!lmdb_store_lookup_qname(st->lmdb, qname_hash(qname), qname, &lr))
        return 0;
    fill_record_from_lmdb(&lr, out);
    lmdb_store_record_clear(&lr);
    return 1;
}

int store_lookup_read_id(store_t *st, int64_t read_id, store_record_t *out) {
    lmdb_store_record_t lr;
    if (!lmdb_store_lookup_read_id(st->lmdb, read_id, &lr))
        return 0;
    fill_record_from_lmdb(&lr, out);
    lmdb_store_record_clear(&lr);
    return 1;
}

void store_record_clear(store_record_t *rec) {
    if (!rec) return;
    free(rec->qname);
    memset(rec, 0, sizeof *rec);
}

/* Fetch payload by BGZF voffset (BAM only). */
static int fetch_bam_voffset(store_t *st, const char *qname,
                             int64_t voffset, bam1_t *out) {
    BGZF *bgzf_fp = st->sam->fp.bgzf;
    if (bgzf_seek(bgzf_fp, voffset, SEEK_SET) < 0)
        die("bgzf_seek to voffset %lld failed in '%s'",
            (long long)voffset, st->store_path);
    int rc = sam_read1(st->sam, st->hdr, out);
    if (rc < 0) return ISAM_ERR_IO;
    const char *got = bam_get_qname(out);
    if (!got || strcmp(got, qname) != 0)
        die("voffset hint stale: expected qname '%s', got '%s' "
            "(rebuild store index with `isamtools index`)",
            qname, got ? got : "<null>");
    return 0;
}

/* Fetch payload by region iteration; verify QNAME. */
static int fetch_iter(store_t *st, const char *qname,
                      int tid, int64_t pos0, bam1_t *out) {
    /* Need the BAM/CRAM index for region queries. */
    static __thread hts_idx_t *cached_idx = NULL;
    static __thread store_t   *cached_st  = NULL;
    if (cached_st != st) {
        if (cached_idx) hts_idx_destroy(cached_idx);
        cached_idx = sam_index_load(st->sam, st->store_path);
        if (!cached_idx)
            die("cannot load index for canonical store '%s' "
                "(samtools index it first)", st->store_path);
        cached_st = st;
    }
    if (tid < 0 || pos0 < 0) {
        /* unmapped: iterate unmapped section */
        hts_itr_t *it = sam_itr_queryi(cached_idx, HTS_IDX_NOCOOR, 0, 0);
        if (!it) return ISAM_ERR_IO;
        int rc;
        while ((rc = sam_itr_next(st->sam, it, out)) >= 0) {
            if (strcmp(bam_get_qname(out), qname) == 0) {
                hts_itr_destroy(it);
                return 0;
            }
        }
        hts_itr_destroy(it);
        return ISAM_ERR_MISSING_QNAME;
    }
    hts_itr_t *it = sam_itr_queryi(cached_idx, tid, pos0, pos0 + 1);
    if (!it) return ISAM_ERR_IO;
    int rc;
    while ((rc = sam_itr_next(st->sam, it, out)) >= 0) {
        if (out->core.pos == pos0 &&
            strcmp(bam_get_qname(out), qname) == 0) {
            hts_itr_destroy(it);
            return 0;
        }
        if (out->core.pos > pos0) break;
    }
    hts_itr_destroy(it);
    return ISAM_ERR_MISSING_QNAME;
}

int store_fetch_payload(store_t *st, const char *qname,
                        const store_record_t *rec, bam1_t *out) {
    if (st->is_bam && rec->has_voffset && rec->bam_voffset >= 0)
        return fetch_bam_voffset(st, qname, rec->bam_voffset, out);
    return fetch_iter(st, qname,
                      rec->is_unmapped ? -1 : rec->tid,
                      rec->is_unmapped ? -1 : rec->pos0,
                      out);
}

payload_cache_t *payload_cache_create(size_t max_bytes) {
    payload_cache_t *cache = xcalloc(1, sizeof *cache);
    cache->max_bytes = max_bytes ? max_bytes : 64 * 1024 * 1024;
    cache->cap_entries = cache->max_bytes / 4096;
    if (cache->cap_entries < 16) cache->cap_entries = 16;
    if (cache->cap_entries > 4096) cache->cap_entries = 4096;
    cache->entries = xcalloc(cache->cap_entries, sizeof *cache->entries);
    size_t buckets = 1;
    while (buckets < cache->cap_entries * 2) buckets <<= 1;
    cache->bucket_mask = buckets - 1;
    cache->buckets = xmalloc(buckets * sizeof *cache->buckets);
    for (size_t i = 0; i < buckets; ++i) cache->buckets[i] = -1;
    return cache;
}

void payload_cache_destroy(payload_cache_t *cache) {
    if (!cache) return;
    for (size_t i = 0; i < cache->n_entries; ++i)
        if (cache->entries[i].bam) bam_destroy1(cache->entries[i].bam);
    free(cache->buckets);
    free(cache->entries);
    free(cache);
}

static size_t payload_cache_bucket(payload_cache_t *cache, int64_t read_id) {
    uint64_t value = (uint64_t)read_id;
    value ^= value >> 33;
    value *= 0xff51afd7ed558ccdULL;
    value ^= value >> 33;
    return (size_t)value & cache->bucket_mask;
}

static int payload_cache_find(payload_cache_t *cache, int64_t read_id) {
    size_t bucket = payload_cache_bucket(cache, read_id);
    for (int idx = cache->buckets[bucket]; idx >= 0;
         idx = cache->entries[idx].next) {
        if (cache->entries[idx].read_id == read_id)
            return idx;
    }
    return -1;
}

static void payload_cache_link(payload_cache_t *cache, size_t slot) {
    size_t bucket = payload_cache_bucket(cache, cache->entries[slot].read_id);
    cache->entries[slot].next = cache->buckets[bucket];
    cache->buckets[bucket] = (int)slot;
}

static void payload_cache_unlink(payload_cache_t *cache, size_t slot) {
    if (!cache->entries[slot].bam) return;
    size_t bucket = payload_cache_bucket(cache, cache->entries[slot].read_id);
    int *cur = &cache->buckets[bucket];
    while (*cur >= 0) {
        if ((size_t)*cur == slot) {
            *cur = cache->entries[*cur].next;
            return;
        }
        cur = &cache->entries[*cur].next;
    }
}

static void payload_cache_evict(payload_cache_t *cache, size_t slot) {
    if (!cache->entries[slot].bam) return;
    payload_cache_unlink(cache, slot);
    cache->used_bytes -= bam_payload_bytes(cache->entries[slot].bam);
    bam_destroy1(cache->entries[slot].bam);
    cache->entries[slot].bam = NULL;
    cache->entries[slot].next = -1;
}

static size_t payload_cache_slot(payload_cache_t *cache) {
    if (cache->n_entries < cache->cap_entries)
        return cache->n_entries++;

    for (size_t tries = 0; tries < cache->cap_entries; ++tries) {
        size_t slot = cache->evict_cursor++ % cache->cap_entries;
        payload_cache_evict(cache, slot);
        return slot;
    }
    return 0;
}

static void payload_cache_insert(payload_cache_t *cache,
                                 int64_t read_id, const bam1_t *b) {
    size_t bytes = bam_payload_bytes(b);
    if (bytes > cache->max_bytes) return;

    size_t slot = payload_cache_slot(cache);

    cache->entries[slot].read_id = read_id;
    cache->entries[slot].bam = bam_dup1(b);
    if (!cache->entries[slot].bam) die("bam_dup1 failed");
    cache->used_bytes += bam_payload_bytes(cache->entries[slot].bam);
    payload_cache_link(cache, slot);

    for (size_t tries = 0;
         cache->used_bytes > cache->max_bytes && tries < cache->cap_entries;
         ++tries) {
        size_t victim = cache->evict_cursor++ % cache->cap_entries;
        if (victim != slot) payload_cache_evict(cache, victim);
    }
}

int payload_cache_fetch(payload_cache_t *cache, store_t *st,
                        const store_record_t *rec, bam1_t *out) {
    int idx = payload_cache_find(cache, rec->read_id);
    if (idx >= 0) {
        if (!bam_copy1(out, cache->entries[idx].bam))
            die("bam_copy1 failed");
        return 0;
    }

    int rc = store_fetch_payload(st, rec->qname, rec, out);
    if (rc != 0) return rc;
    payload_cache_insert(cache, rec->read_id, out);
    return 0;
}
