#include "store.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/cram.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>

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

/* ---- fetch instrumentation (enable with ISAM_FETCH_STATS=1) ---------- */
static long long g_fs_voffset = 0;
static long long g_fs_cramloc_call = 0, g_fs_cramloc_hit = 0, g_fs_cramloc_decode = 0;
static long long g_fs_cram_decode_recs = 0;
static long long g_fs_iter_call = 0, g_fs_iter_mapped = 0, g_fs_iter_unmapped = 0;

/* ---- CRAM decoded-slice LRU cache -------------------------------------
 * CRAM random access is slice/container granular: fetching one read means
 * seeking to its container and decoding the whole slice. This cache holds a
 * bounded set of decoded containers keyed by byte offset, so consecutive reads
 * from the same container decode it once, not once per read. Eviction is
 * least-recently-used down to a byte budget. */

typedef struct cram_slice {
    int64_t   coffset;    /* container byte offset; -1 = free slot */
    bam1_t  **recs;       /* recs[0..nrec-1], owned */
    int       nrec;
    size_t    bytes;
    uint64_t  last_used;
} cram_slice_t;

struct cram_slice_cache {
    cram_slice_t *slots;
    int      cap;
    size_t   max_bytes;
    size_t   used_bytes;
    uint64_t clock;
};

static cram_slice_cache_t *cram_slice_cache_create(size_t max_bytes) {
    cram_slice_cache_t *c = xcalloc(1, sizeof *c);
    c->max_bytes = max_bytes ? max_bytes : (size_t)256 * 1024 * 1024;
    c->cap = 4096;
    c->slots = xcalloc((size_t)c->cap, sizeof *c->slots);
    for (int i = 0; i < c->cap; ++i) c->slots[i].coffset = -1;
    return c;
}

static void cram_slice_free(cram_slice_t *s) {
    if (s->recs) {
        for (int i = 0; i < s->nrec; ++i)
            if (s->recs[i]) bam_destroy1(s->recs[i]);
        free(s->recs);
        s->recs = NULL;
    }
    s->nrec = 0;
    s->coffset = -1;
    s->bytes = 0;
}

static void cram_slice_cache_destroy(cram_slice_cache_t *c) {
    if (!c) return;
    for (int i = 0; i < c->cap; ++i) cram_slice_free(&c->slots[i]);
    free(c->slots);
    free(c);
}

static int cram_slice_cache_find(cram_slice_cache_t *c, int64_t coffset) {
    for (int i = 0; i < c->cap; ++i)
        if (c->slots[i].coffset == coffset) return i;
    return -1;
}

/* A free slot, else evict the least-recently-used occupied slot and reuse it. */
static int cram_slice_cache_victim(cram_slice_cache_t *c) {
    int lru = -1;
    uint64_t best = UINT64_MAX;
    for (int i = 0; i < c->cap; ++i) {
        if (c->slots[i].coffset < 0) return i;
        if (c->slots[i].last_used < best) { best = c->slots[i].last_used; lru = i; }
    }
    c->used_bytes -= c->slots[lru].bytes;
    cram_slice_free(&c->slots[lru]);
    return lru;
}

/* nrec for a container, from the (ascending) offset table. */
static int cram_container_nrec(store_t *st, int64_t coffset) {
    int64_t lo = 0, hi = st->cram_nc - 1;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (st->cram_coff[mid] == coffset) return st->cram_nrec[mid];
        if (st->cram_coff[mid] < coffset) lo = mid + 1; else hi = mid - 1;
    }
    die("cram container offset %lld not in table for '%s'",
        (long long)coffset, st->store_path);
    return 0;
}

static void cram_slice_decode(store_t *st, int64_t coffset, int nrec,
                              cram_slice_t *slot) {
    cram_fd *fd = st->sam->fp.cram;
    if (cram_seek(fd, (off_t)coffset, SEEK_SET) < 0)
        die("cram_seek to container %lld failed in '%s'",
            (long long)coffset, st->store_path);
    slot->recs = xmalloc((size_t)nrec * sizeof *slot->recs);
    slot->nrec = nrec;
    slot->coffset = coffset;
    slot->bytes = 0;
    for (int i = 0; i < nrec; ++i) {
        bam1_t *r = bam_init1();
        int rc = sam_read1(st->sam, st->hdr, r);
        if (rc < 0)
            die("sam_read1 (%d/%d) in CRAM container %lld of '%s' failed (rc=%d); "
                "is the reference set (-T)?", i, nrec, (long long)coffset,
                st->store_path, rc);
        slot->recs[i] = r;
        slot->bytes += sizeof(bam1_t) + (size_t)r->m_data;
    }
}

/* Fetch payload by CRAM container locator (container offset + within index). */
static int fetch_cram_loc(store_t *st, const char *qname,
                          int64_t coffset, int32_t within, bam1_t *out) {
    cram_slice_cache_t *c = st->slice_cache;
    ++g_fs_cramloc_call;
    int idx = cram_slice_cache_find(c, coffset);
    if (idx >= 0) ++g_fs_cramloc_hit;
    if (idx < 0) {
        ++g_fs_cramloc_decode;
        int nrec = cram_container_nrec(st, coffset);
        g_fs_cram_decode_recs += nrec;
        if (within < 0 || within >= nrec)
            die("cram within-index %d out of range [0,%d) for container %lld in '%s'",
                within, nrec, (long long)coffset, st->store_path);
        idx = cram_slice_cache_victim(c);
        cram_slice_decode(st, coffset, nrec, &c->slots[idx]);
        c->used_bytes += c->slots[idx].bytes;
        while (c->used_bytes > c->max_bytes) {
            int v = -1;
            uint64_t best = UINT64_MAX;
            for (int i = 0; i < c->cap; ++i) {
                if (i == idx || c->slots[i].coffset < 0) continue;
                if (c->slots[i].last_used < best) { best = c->slots[i].last_used; v = i; }
            }
            if (v < 0) break;
            c->used_bytes -= c->slots[v].bytes;
            cram_slice_free(&c->slots[v]);
        }
    }
    cram_slice_t *s = &c->slots[idx];
    s->last_used = ++c->clock;
    if (within < 0 || within >= s->nrec)
        die("cram within-index %d >= %d records in container %lld of '%s'",
            within, s->nrec, (long long)coffset, st->store_path);
    bam1_t *r = s->recs[within];
    const char *got = bam_get_qname(r);
    if (!got || strcmp(got, qname) != 0)
        die("cram locator stale: expected qname '%s', got '%s' in '%s' "
            "(rebuild store index with `isamtools index`)",
            qname, got ? got : "<null>", st->store_path);
    if (!bam_copy1(out, r)) die("bam_copy1 failed");
    return 0;
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

    /* CRAM canonical: load the container table and set up the slice cache.
     * A store indexed before v2 has no table; such CRAM stores fall back to
     * the coordinate-iteration payload fetch. */
    if (!st->is_bam &&
        lmdb_store_get_cram_container_table(st->lmdb, &st->cram_coff,
                                            &st->cram_nrec, &st->cram_nc)) {
        /* Decoded-slice cache budget. A PacBio CRAM container decodes to
         * ~250 MB (10k reads x SEQ/QUAL/methylation), so the old 256 MB default
         * held ~1 container and thrashed: supplementary/secondary alignments
         * verify against the canonical payload stored at the *primary's*
         * coordinate, which back-jumps to distant containers and forced a full
         * re-decode per read (~40x read amplification observed). Default to a
         * large budget (nodes here have 100s of GB RAM); tune with
         * ISAM_SLICE_CACHE_MB. */
        size_t cache_mb = 8192;
        const char *env = getenv("ISAM_SLICE_CACHE_MB");
        if (env && *env) {
            long v = strtol(env, NULL, 10);
            if (v > 0) cache_mb = (size_t)v;
        }
        st->slice_cache = cram_slice_cache_create(cache_mb * 1024 * 1024);
    }

    return st;
}

void store_close(store_t *st) {
    if (!st) return;
    if (getenv("ISAM_FETCH_STATS")) {
        long long cl = g_fs_cramloc_call;
        fprintf(stderr,
            "isamtools fetch-stats: voffset=%lld cram_loc=%lld (hit=%lld miss/decode=%lld, "
            "hit%%=%.1f, decoded_recs=%lld) iter=%lld (mapped=%lld unmapped=%lld)\n",
            g_fs_voffset, cl, g_fs_cramloc_hit, g_fs_cramloc_decode,
            cl ? 100.0 * (double)g_fs_cramloc_hit / (double)cl : 0.0,
            g_fs_cram_decode_recs, g_fs_iter_call, g_fs_iter_mapped, g_fs_iter_unmapped);
    }
    cram_slice_cache_destroy(st->slice_cache);
    free(st->cram_coff);
    free(st->cram_nrec);
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
    out->cram_coffset = lr->cram_coffset;
    out->cram_within = lr->cram_within;
    out->has_cram_loc = lr->has_cram_loc;
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
    ++g_fs_iter_call;
    if (tid < 0 || pos0 < 0) ++g_fs_iter_unmapped; else ++g_fs_iter_mapped;
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
    if (st->is_bam && rec->has_voffset && rec->bam_voffset >= 0) {
        ++g_fs_voffset;
        return fetch_bam_voffset(st, qname, rec->bam_voffset, out);
    }
    if (!st->is_bam && rec->has_cram_loc && rec->cram_coffset >= 0 && st->slice_cache)
        return fetch_cram_loc(st, qname, rec->cram_coffset, rec->cram_within, out);
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
