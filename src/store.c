#include "store.h"
#include "aux_tags.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/cram.h>
#include <htslib/thread_pool.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <limits.h>
#include <inttypes.h>
#include <pthread.h>
#include <time.h>

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
static double g_fs_cram_decode_seconds = 0;

static double monotonic_seconds(void) {
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) die("clock_gettime failed");
    return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

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
    int       retained;
    size_t    bytes;
    int       prev, next;
    struct cram_decode_task *pending; /* main-thread-owned, outside the LRU */
} cram_slice_t;

struct cram_slice_cache {
    cram_slice_t *slots;
    int      cap;
    size_t   max_bytes;
    size_t   used_bytes;
    size_t   peak_bytes;
    int      head, tail, last;
    uint32_t **remaining_uses;
    uint64_t planned_requests, remaining_requests, released_records;
    size_t   plan_bytes;
    double   plan_started, plan_seconds;
    int      planning;
};

struct cram_decoder_pool {
    hts_tpool *pool;
    hts_tpool_process *queue;
    samFile **files;
    sam_hdr_t **headers;
    struct cram_decode_task **active, **pending;
    pthread_mutex_t mutex;
    pthread_cond_t completed;
    int count;
    int pending_limit, pending_count, peak_pending;
    long long submitted, collected, waits;
    double wait_seconds;
};

static cram_slice_cache_t *cram_slice_cache_create(size_t max_bytes, int64_t containers) {
    if (containers < 1 || containers > INT_MAX)
        die("invalid CRAM container count: %lld", (long long)containers);
    cram_slice_cache_t *c = xcalloc(1, sizeof *c);
    c->max_bytes = max_bytes ? max_bytes : (size_t)256 * 1024 * 1024;
    /* One small descriptor per indexed container, but allocate BAM payloads
     * only while resident. Container ordinal is the cache key: no slot scan. */
    c->cap = (int)containers;
    c->head = c->tail = c->last = -1;
    c->slots = xcalloc((size_t)c->cap, sizeof *c->slots);
    for (int i = 0; i < c->cap; ++i) {
        c->slots[i].coffset = -1;
        c->slots[i].prev = c->slots[i].next = -1;
    }
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
    s->retained = 0;
    s->coffset = -1;
    s->bytes = 0;
}

static void cram_slice_cache_destroy(cram_slice_cache_t *c) {
    if (!c) return;
    for (int i = 0; i < c->cap; ++i) {
        cram_slice_free(&c->slots[i]);
        if (c->remaining_uses) free(c->remaining_uses[i]);
    }
    free(c->remaining_uses);
    free(c->slots);
    free(c);
}

static void cram_slice_unlink(cram_slice_cache_t *c, int idx) {
    cram_slice_t *s = &c->slots[idx];
    if (s->prev >= 0) c->slots[s->prev].next = s->next;
    else c->head = s->next;
    if (s->next >= 0) c->slots[s->next].prev = s->prev;
    else c->tail = s->prev;
    s->prev = s->next = -1;
}

static void cram_slice_touch(cram_slice_cache_t *c, int idx, int resident) {
    if (c->tail == idx) return;
    if (resident) cram_slice_unlink(c, idx);
    cram_slice_t *s = &c->slots[idx];
    s->prev = c->tail;
    if (c->tail >= 0) c->slots[c->tail].next = idx;
    else c->head = idx;
    c->tail = idx;
}

/* Container ordinal from the immutable, ascending locator table. */
static int cram_container_index(store_t *st, int64_t coffset) {
    cram_slice_cache_t *c = st->slice_cache;
    if (c->last >= 0 && st->cram_coff[c->last] == coffset) return c->last;
    int64_t lo = 0, hi = st->cram_nc - 1;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (st->cram_coff[mid] == coffset) return c->last = (int)mid;
        if (st->cram_coff[mid] < coffset) lo = mid + 1; else hi = mid - 1;
    }
    die("cram container offset %lld not in table for '%s'",
        (long long)coffset, st->store_path);
    return 0;
}

static void cram_slice_decode(samFile *input, sam_hdr_t *header, const char *path,
                              int64_t coffset, int nrec, cram_slice_t *slot,
                              const uint32_t *remaining_uses) {
    cram_fd *fd = input->fp.cram;
    if (cram_seek(fd, (off_t)coffset, SEEK_SET) < 0)
        die("cram_seek to container %lld failed in '%s'",
            (long long)coffset, path);
    slot->recs = xcalloc((size_t)nrec, sizeof *slot->recs);
    slot->nrec = nrec;
    slot->coffset = coffset;
    slot->bytes = 0;
    slot->retained = 0;
    bam1_t *discard = remaining_uses ? bam_init1() : NULL;
    for (int i = 0; i < nrec; ++i) {
        int keep = !remaining_uses || remaining_uses[i] != 0;
        bam1_t *r = keep ? bam_init1() : discard;
        int rc = sam_read1(input, header, r);
        if (rc < 0)
            die("sam_read1 (%d/%d) in CRAM container %lld of '%s' failed (rc=%d); "
                "is the reference set (-T)?", i, nrec, (long long)coffset,
                path, rc);
        if (keep) {
            slot->recs[i] = r;
            slot->bytes += bam_payload_bytes(r);
            ++slot->retained;
        }
    }
    bam_destroy1(discard);
}

void store_payload_plan_begin(store_t *st) {
    cram_slice_cache_t *c = st->slice_cache;
    if (!c || c->remaining_uses || c->used_bytes ||
        (st->decoders && st->decoders->submitted))
        die("payload planning requires an unused indexed CRAM cache");
    c->remaining_uses = xcalloc((size_t)c->cap, sizeof *c->remaining_uses);
    c->plan_bytes = (size_t)c->cap * sizeof *c->remaining_uses;
    c->planning = 1;
    c->plan_started = monotonic_seconds();
}

void store_payload_plan_add(store_t *st, const store_record_t *rec) {
    cram_slice_cache_t *c = st->slice_cache;
    if (!c || !c->planning || !rec->has_cram_loc)
        die("payload planning requires a CRAM locator for every selected record");
    int idx = cram_container_index(st, rec->cram_coffset);
    int nrec = st->cram_nrec[idx];
    if (rec->cram_within < 0 || rec->cram_within >= nrec)
        die("invalid CRAM within-index in payload plan");
    if (!c->remaining_uses[idx]) {
        c->remaining_uses[idx] = xcalloc((size_t)nrec, sizeof **c->remaining_uses);
        c->plan_bytes += (size_t)nrec * sizeof **c->remaining_uses;
    }
    uint32_t *uses = &c->remaining_uses[idx][rec->cram_within];
    if (*uses == UINT32_MAX) die("payload plan request count overflow");
    ++*uses;
    ++c->planned_requests;
    ++c->remaining_requests;
}

void store_payload_plan_finish(store_t *st) {
    cram_slice_cache_t *c = st->slice_cache;
    if (!c || !c->planning) die("no payload plan is being built");
    c->planning = 0;
    c->plan_seconds = monotonic_seconds() - c->plan_started;
}

void store_payload_release(store_t *st, const store_record_t *rec) {
    cram_slice_cache_t *c = st->slice_cache;
    if (!c || !c->remaining_uses) return;
    int idx = cram_container_index(st, rec->cram_coffset);
    cram_slice_t *slot = &c->slots[idx];
    if (c->planning || !rec->has_cram_loc || rec->cram_within < 0 ||
        rec->cram_within >= slot->nrec || !slot->recs ||
        !c->remaining_uses[idx] || !c->remaining_uses[idx][rec->cram_within] ||
        !slot->recs[rec->cram_within])
        die("payload release does not match the coordinate request plan");
    --c->remaining_requests;
    if (--c->remaining_uses[idx][rec->cram_within]) return;
    bam1_t *payload = slot->recs[rec->cram_within];
    size_t bytes = bam_payload_bytes(payload);
    slot->bytes -= bytes;
    c->used_bytes -= bytes;
    bam_destroy1(payload);
    slot->recs[rec->cram_within] = NULL;
    ++c->released_records;
    if (--slot->retained == 0) {
        cram_slice_unlink(c, idx);
        cram_slice_free(slot);
    }
}

static void cram_slice_install(cram_slice_cache_t *c, int idx) {
    c->used_bytes += c->slots[idx].bytes;
    if (c->used_bytes > c->peak_bytes) c->peak_bytes = c->used_bytes;
    cram_slice_touch(c, idx, 0);
    /* Keep the requested container even if it alone exceeds the budget. */
    while (c->used_bytes > c->max_bytes && c->head != idx) {
        int victim = c->head;
        if (victim < 0) die("invalid CRAM cache LRU state");
        cram_slice_unlink(c, victim);
        c->used_bytes -= c->slots[victim].bytes;
        cram_slice_free(&c->slots[victim]);
    }
}

typedef struct cram_decode_task {
    samFile *input;
    sam_hdr_t *header;
    const char *path;
    int idx, nrec;
    int64_t coffset;
    cram_slice_t decoded;
    double seconds;
    const uint32_t *remaining_uses;
    cram_decoder_pool_t *decoders;
    int ticket, done, ready;
} cram_decode_task_t;

static void *cram_decode_task(void *argument) {
    cram_decode_task_t *task = argument;
    double started = monotonic_seconds();
    cram_slice_decode(task->input, task->header, task->path, task->coffset,
                      task->nrec, &task->decoded, task->remaining_uses);
    task->seconds = monotonic_seconds() - started;
    cram_decoder_pool_t *decoders = task->decoders;
    if (pthread_mutex_lock(&decoders->mutex) != 0) die("cannot lock CRAM decoder pool");
    task->done = 1;
    if (pthread_cond_broadcast(&decoders->completed) != 0) die("cannot signal CRAM decode completion");
    if (pthread_mutex_unlock(&decoders->mutex) != 0) die("cannot unlock CRAM decoder pool");
    return NULL;
}

/* Only the main thread collects results or changes cache residency. Waiting
 * is for one demanded container, never all the jobs in a prefetch batch. */
static void cram_decoder_collect(cram_decoder_pool_t *decoders, cram_decode_task_t *needed) {
    if (pthread_mutex_lock(&decoders->mutex) != 0) die("cannot lock CRAM decoder pool");
    if (needed && !needed->done) {
        double started = monotonic_seconds();
        ++decoders->waits;
        while (!needed->done)
            if (pthread_cond_wait(&decoders->completed, &decoders->mutex) != 0)
                die("cannot wait for CRAM container decode");
        decoders->wait_seconds += monotonic_seconds() - started;
    }
    for (int i = 0; i < decoders->count; ++i) {
        cram_decode_task_t *task = decoders->active[i];
        if (!task || !task->done) continue;
        task->ready = 1;
        decoders->active[i] = NULL;
        ++decoders->collected;
        ++g_fs_cramloc_decode;
        g_fs_cram_decode_recs += task->nrec;
        g_fs_cram_decode_seconds += task->seconds;
    }
    if (pthread_mutex_unlock(&decoders->mutex) != 0) die("cannot unlock CRAM decoder pool");
}

void store_set_decode_threads(store_t *st, int threads, const char *reference) {
    if (threads < 1 || threads > 64 || st->is_bam || !st->slice_cache || st->decoders)
        die("container decode threads require an indexed CRAM and 1..64 workers");
    cram_decoder_pool_t *decoders = xcalloc(1, sizeof *decoders);
    decoders->count = threads;
    decoders->pending_limit = threads * 2;
    decoders->active = xcalloc((size_t)threads, sizeof *decoders->active);
    decoders->pending = xcalloc((size_t)decoders->pending_limit, sizeof *decoders->pending);
    if (pthread_mutex_init(&decoders->mutex, NULL) != 0 ||
        pthread_cond_init(&decoders->completed, NULL) != 0)
        die("cannot initialize CRAM decoder synchronization");
    decoders->files = xcalloc((size_t)threads, sizeof *decoders->files);
    decoders->headers = xcalloc((size_t)threads, sizeof *decoders->headers);
    for (int i = 0; i < threads; ++i) {
        decoders->files[i] = sam_open(st->store_path, "r");
        if (!decoders->files[i]) die("cannot open CRAM decoder '%s'", st->store_path);
        decoders->headers[i] = sam_hdr_read(decoders->files[i]);
        if (!decoders->headers[i]) die("cannot read CRAM decoder header '%s'", st->store_path);
        if (reference && hts_set_fai_filename(decoders->files[i], reference) != 0)
            die("cannot set CRAM decoder reference '%s'", reference);
    }
    decoders->pool = hts_tpool_init(threads);
    if (!decoders->pool) die("cannot create CRAM container decoder pool");
    /* Completion belongs to each task rather than HTSlib's ordered result
     * queue, so a slow speculative task cannot block a demanded result. */
    decoders->queue = hts_tpool_process_init(decoders->pool, threads, 1);
    if (!decoders->queue) die("cannot create CRAM container decoder queue");
    st->decoders = decoders;
}

int store_payload_cached(store_t *st, const store_record_t *rec) {
    if (!rec->has_cram_loc || !st->slice_cache) return 0;
    return st->slice_cache->slots[cram_container_index(st, rec->cram_coffset)].recs != NULL;
}

size_t store_prefetch_payloads(store_t *st, const store_record_t *records, size_t count) {
    cram_decoder_pool_t *decoders = st->decoders;
    if (!decoders) return count;
    cram_decoder_collect(decoders, NULL);
    for (size_t i = 0; i < count; ++i) {
        if (!records[i].has_cram_loc) continue;
        int idx = cram_container_index(st, records[i].cram_coffset);
        cram_slice_t *slot = &st->slice_cache->slots[idx];
        if (slot->recs || slot->pending) continue;
        if (decoders->pending_count == decoders->pending_limit) return i;
        int worker = 0, ticket = 0;
        while (worker < decoders->count && decoders->active[worker]) ++worker;
        if (worker == decoders->count) return i;
        while (ticket < decoders->pending_limit && decoders->pending[ticket]) ++ticket;
        if (ticket == decoders->pending_limit) die("invalid CRAM prefetch ticket state");
        cram_decode_task_t *task = xcalloc(1, sizeof *task);
        task->input = decoders->files[worker];
        task->header = decoders->headers[worker];
        task->path = st->store_path;
        task->idx = idx;
        task->ticket = ticket;
        task->decoders = decoders;
        task->coffset = st->cram_coff[idx];
        task->nrec = st->cram_nrec[idx];
        if (st->slice_cache->remaining_uses) {
            if (st->slice_cache->planning || !st->slice_cache->remaining_uses[idx])
                die("prefetch is missing a completed payload plan");
            task->remaining_uses = st->slice_cache->remaining_uses[idx];
        }
        task->decoded.prev = task->decoded.next = -1;
        slot->pending = task;
        decoders->active[worker] = decoders->pending[ticket] = task;
        ++decoders->submitted;
        if (++decoders->pending_count > decoders->peak_pending)
            decoders->peak_pending = decoders->pending_count;
        if (hts_tpool_dispatch(decoders->pool, decoders->queue, cram_decode_task, task) != 0)
            die("cannot submit CRAM container decode task");
    }
    return count;
}

/* Fetch payload by CRAM container locator (container offset + within index). */
static const bam1_t *borrow_cram_loc(store_t *st, const char *qname,
                                    int64_t coffset, int32_t within) {
    cram_slice_cache_t *c = st->slice_cache;
    ++g_fs_cramloc_call;
    int idx = cram_container_index(st, coffset);
    int resident = c->slots[idx].recs != NULL;
    if (resident) {
        ++g_fs_cramloc_hit;
        cram_slice_touch(c, idx, 1);
    } else if (c->slots[idx].pending) {
        cram_decode_task_t *task = c->slots[idx].pending;
        cram_decoder_collect(st->decoders, task);
        if (!task->ready) die("missing demanded CRAM decode result");
        c->slots[idx] = task->decoded;
        st->decoders->pending[task->ticket] = NULL;
        --st->decoders->pending_count;
        free(task);
        cram_slice_install(c, idx);
    } else {
        ++g_fs_cramloc_decode;
        int nrec = st->cram_nrec[idx];
        g_fs_cram_decode_recs += nrec;
        if (within < 0 || within >= nrec)
            die("cram within-index %d out of range [0,%d) for container %lld in '%s'",
                within, nrec, (long long)coffset, st->store_path);
        double started = monotonic_seconds();
        cram_slice_decode(st->sam, st->hdr, st->store_path, coffset, nrec, &c->slots[idx],
                          c->remaining_uses ? c->remaining_uses[idx] : NULL);
        g_fs_cram_decode_seconds += monotonic_seconds() - started;
        cram_slice_install(c, idx);
    }
    cram_slice_t *s = &c->slots[idx];
    if (within < 0 || within >= s->nrec)
        die("cram within-index %d >= %d records in container %lld of '%s'",
            within, s->nrec, (long long)coffset, st->store_path);
    bam1_t *r = s->recs[within];
    if (!r) die("requested CRAM payload is absent from the coordinate plan");
    const char *got = bam_get_qname(r);
    if (!got || strcmp(got, qname) != 0)
        die("cram locator stale: expected qname '%s', got '%s' in '%s' "
            "(rebuild store index with `isamtools index`)",
            qname, got ? got : "<null>", st->store_path);
    return r;
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
        st->slice_cache = cram_slice_cache_create(cache_mb * 1024 * 1024, st->cram_nc);
    }

    return st;
}

void store_close(store_t *st) {
    if (!st) return;
    if (st->decoders) {
        if (hts_tpool_process_flush(st->decoders->queue) != 0)
            die("cannot drain CRAM decoder queue");
        cram_decoder_collect(st->decoders, NULL);
    }
    if (getenv("ISAM_FETCH_STATS")) {
        long long cl = g_fs_cramloc_call;
        fprintf(stderr,
            "isamtools fetch-stats: voffset=%lld cram_loc=%lld (hit=%lld miss/decode=%lld, "
            "hit%%=%.1f, decoded_recs=%lld) iter=%lld (mapped=%lld unmapped=%lld)\n",
            g_fs_voffset, cl, g_fs_cramloc_hit, g_fs_cramloc_decode,
            cl ? 100.0 * (double)g_fs_cramloc_hit / (double)cl : 0.0,
            g_fs_cram_decode_recs, g_fs_iter_call, g_fs_iter_mapped, g_fs_iter_unmapped);
        fprintf(stderr, "isamtools decode-stats: seconds=%.6f cache_peak_bytes=%zu\n",
                g_fs_cram_decode_seconds, st->slice_cache ? st->slice_cache->peak_bytes : 0);
        cram_slice_cache_t *c = st->slice_cache;
        if (c && c->remaining_uses)
            fprintf(stderr, "isamtools selection-stats: requests=%" PRIu64
                    " remaining=%" PRIu64 " released_records=%" PRIu64
                    " plan_bytes=%zu plan_seconds=%.6f\n",
                    c->planned_requests, c->remaining_requests, c->released_records,
                    c->plan_bytes, c->plan_seconds);
    }
    if (st->decoders) {
        cram_decoder_pool_t *decoders = st->decoders;
        if (getenv("ISAM_FETCH_STATS"))
            fprintf(stderr, "isamtools pipeline-stats: submitted=%lld completed=%lld "
                    "waits=%lld wait_seconds=%.6f peak_pending=%d pending_limit=%d\n",
                    decoders->submitted, decoders->collected, decoders->waits,
                    decoders->wait_seconds, decoders->peak_pending, decoders->pending_limit);
        hts_tpool_process_destroy(decoders->queue);
        hts_tpool_destroy(decoders->pool);
        for (int i = 0; i < decoders->pending_limit; ++i) {
            cram_decode_task_t *task = decoders->pending[i];
            if (!task) continue;
            st->slice_cache->slots[task->idx].pending = NULL;
            cram_slice_free(&task->decoded);
            free(task);
        }
        for (int i = 0; i < decoders->count; ++i) {
            sam_hdr_destroy(decoders->headers[i]);
            if (sam_close(decoders->files[i]) != 0) die("cannot close CRAM decoder");
        }
        free(decoders->headers);
        free(decoders->files);
        free(decoders->active);
        free(decoders->pending);
        if (pthread_cond_destroy(&decoders->completed) != 0 ||
            pthread_mutex_destroy(&decoders->mutex) != 0)
            die("cannot destroy CRAM decoder synchronization");
        free(decoders);
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
    out->has_tag_fingerprints = lr->has_tag_fingerprints;
    out->n_tag_fingerprints = lr->n_tag_fingerprints;
    if (lr->n_tag_fingerprints > 0) {
        out->tag_fingerprints = xmalloc(lr->n_tag_fingerprints *
                                        sizeof *out->tag_fingerprints);
        memcpy(out->tag_fingerprints, lr->tag_fingerprints,
               lr->n_tag_fingerprints * sizeof *out->tag_fingerprints);
    }
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
    free(rec->tag_fingerprints);
    memset(rec, 0, sizeof *rec);
}

int store_record_tag_fingerprint(const store_record_t *rec,
                                 const char tag[2],
                                 uint64_t *fingerprint) {
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < rec->n_tag_fingerprints; ++i) {
        if (rec->tag_fingerprints[i].tag_code != code)
            continue;
        *fingerprint = rec->tag_fingerprints[i].fingerprint;
        return 1;
    }
    return 0;
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
    if (!st->is_bam && rec->has_cram_loc && rec->cram_coffset >= 0 && st->slice_cache) {
        const bam1_t *payload = borrow_cram_loc(st, qname, rec->cram_coffset, rec->cram_within);
        if (!bam_copy1(out, payload)) die("bam_copy1 failed");
        return 0;
    }
    return fetch_iter(st, qname,
                      rec->is_unmapped ? -1 : rec->tid,
                      rec->is_unmapped ? -1 : rec->pos0,
                      out);
}

int store_borrow_payload(store_t *st, const store_record_t *rec,
                         bam1_t *scratch, const bam1_t **payload) {
    if (!st->is_bam && rec->has_cram_loc && rec->cram_coffset >= 0 && st->slice_cache) {
        *payload = borrow_cram_loc(st, rec->qname, rec->cram_coffset, rec->cram_within);
        return 0;
    }
    int rc = store_fetch_payload(st, rec->qname, rec, scratch);
    *payload = rc == 0 ? scratch : NULL;
    return rc;
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
