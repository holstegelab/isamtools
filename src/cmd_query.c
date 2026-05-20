#include "common.h"
#include "backend_lmdb.h"
#include "store.h"
#include "view.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <errno.h>
#include <libgen.h>
#include <sys/stat.h>

typedef enum query_source_type {
    QUERY_SOURCE_STORE = 1,
    QUERY_SOURCE_VIEW = 2
} query_source_type_t;

typedef struct qname_set {
    char **items;
    size_t n;
    size_t cap;
    int sorted;
} qname_set_t;

typedef struct query_aln {
    int tid;
    int64_t pos0;
    int64_t end0;
    int flag;
    int mapq;
    char *cigar;
    int is_unmapped;
} query_aln_t;

typedef struct query_hits {
    int has_read;
    int64_t read_id;
    query_aln_t *alns;
    size_t n_alns;
    size_t cap_alns;
} query_hits_t;

typedef struct query_interval {
    int tid;
    int64_t beg;
    int64_t end;
} query_interval_t;

typedef struct query_source {
    char *alias;
    char *path;
    query_source_type_t type;
    store_t *store;
    view_t *view;
    sam_hdr_t *hdr;
    int needs_payload;
} query_source_t;

typedef enum query_column_kind {
    QUERY_COL_QNAME = 1,
    QUERY_COL_SOURCE_FIELD = 2,
    QUERY_COL_SOURCE_SUBSET = 3
} query_column_kind_t;

typedef enum query_field {
    QUERY_FIELD_READ_ID = 1,
    QUERY_FIELD_N = 2,
    QUERY_FIELD_RNAME = 3,
    QUERY_FIELD_POS = 4,
    QUERY_FIELD_START0 = 5,
    QUERY_FIELD_END0 = 6,
    QUERY_FIELD_REGION = 7,
    QUERY_FIELD_FLAG = 8,
    QUERY_FIELD_MAPQ = 9,
    QUERY_FIELD_STRAND = 10,
    QUERY_FIELD_CIGAR = 11,
    QUERY_FIELD_IS_MAPPED = 12
} query_field_t;

typedef struct query_column {
    char *header;
    query_column_kind_t kind;
    int source_index;
    query_field_t field;
    char *subset_name;
} query_column_t;

typedef struct selector_set {
    qname_set_t qnames;
    char *label;
} selector_set_t;

typedef struct span_spec {
    char *arg;
} span_spec_t;

static void usage(void) {
    fputs(
        "Usage: isamtools query [options]\n"
        "\n"
        "Build a read-name selection and emit TSV fields from one or more named\n"
        "store .ssam / alignment-view .isam sources.\n"
        "\n"
        "Sources:\n"
        "  --source A:PATH       Add a named source. PATH may be .ssam, .isam,\n"
        "                        or a store BAM/CRAM/SAM with a matching .ssam.\n"
        "\n"
        "Read selection:\n"
        "  --select A:REGION     Select reads overlapping REGION in source A.\n"
        "  --region A:REGION     Alias for --select.\n"
        "  --subset A:NAME       Select reads in subset NAME from source A.\n"
        "  --qnames FILE         Select one QNAME per line from FILE.\n"
        "  --join union|inner    Combine multiple selectors [union].\n"
        "  --span A:DIST         Expand the selected reads through source A:\n"
        "                        merge selected-read alignments separated by\n"
        "                        <= DIST bases, then add all reads in those spans.\n"
        "\n"
        "Output:\n"
        "  --columns LIST        Comma-separated columns. Examples:\n"
        "                        qname,hg38.region,chm13.region,chm13.in_subset:hifi\n"
        "  --no-header           Do not print the TSV header.\n"
        "  -o, --output FILE     Write TSV to FILE [stdout].\n"
        "  -h, --help            Show this help.\n",
        stderr);
}

static void qname_set_add(qname_set_t *set, const char *qname) {
    if (set->n == set->cap) {
        set->cap = set->cap ? set->cap * 2 : 1024;
        set->items = xrealloc(set->items, set->cap * sizeof *set->items);
    }
    set->items[set->n++] = xstrdup(qname);
    set->sorted = 0;
}

static int qname_ptr_cmp(const void *a, const void *b) {
    const char *const *x = a;
    const char *const *y = b;
    return strcmp(*x, *y);
}

static void qname_set_sort_unique(qname_set_t *set) {
    if (set->sorted)
        return;
    qsort(set->items, set->n, sizeof *set->items, qname_ptr_cmp);
    size_t out = 0;
    for (size_t i = 0; i < set->n; ++i) {
        if (out && strcmp(set->items[out - 1], set->items[i]) == 0) {
            free(set->items[i]);
            continue;
        }
        set->items[out++] = set->items[i];
    }
    set->n = out;
    set->sorted = 1;
}

static int qname_set_contains(qname_set_t *set, const char *qname) {
    qname_set_sort_unique(set);
    return bsearch(&qname, set->items, set->n, sizeof *set->items,
                   qname_ptr_cmp) != NULL;
}

static void qname_set_union_into(qname_set_t *dst, qname_set_t *src) {
    qname_set_sort_unique(src);
    for (size_t i = 0; i < src->n; ++i)
        qname_set_add(dst, src->items[i]);
}

static void qname_set_intersect_into(qname_set_t *dst,
                                     qname_set_t *a,
                                     qname_set_t *b) {
    qname_set_sort_unique(a);
    qname_set_sort_unique(b);
    for (size_t i = 0; i < a->n; ++i)
        if (qname_set_contains(b, a->items[i]))
            qname_set_add(dst, a->items[i]);
}

static void qname_set_destroy(qname_set_t *set) {
    for (size_t i = 0; i < set->n; ++i)
        free(set->items[i]);
    free(set->items);
    memset(set, 0, sizeof *set);
}

static char *xstrndup_(const char *s, size_t n) {
    char *out = xmalloc(n + 1);
    memcpy(out, s, n);
    out[n] = 0;
    return out;
}

static int split_alias_path(const char *arg, char **alias, char **path) {
    const char *colon = strchr(arg, ':');
    if (!colon || colon == arg) {
        *alias = NULL;
        *path = xstrdup(arg);
        return 0;
    }
    const char *slash = strchr(arg, '/');
    if (slash && slash < colon) {
        *alias = NULL;
        *path = xstrdup(arg);
        return 0;
    }
    *alias = xstrndup_(arg, (size_t)(colon - arg));
    *path = xstrdup(colon + 1);
    return 1;
}

static void split_alias_value(const char *arg,
                              const char *what,
                              char **alias,
                              char **value) {
    const char *colon = strchr(arg, ':');
    if (!colon || colon == arg || !colon[1])
        die("%s must be ALIAS:VALUE, got '%s'", what, arg);
    *alias = xstrndup_(arg, (size_t)(colon - arg));
    *value = xstrdup(colon + 1);
}

static char *default_alias_for_path(const char *path, query_source_type_t type) {
    char *copy = xstrdup(path);
    char *base = basename(copy);
    char *alias = xstrdup(base);
    free(copy);
    if (isam_has_path_suffix(alias, ".isam")) {
        alias[strlen(alias) - 5] = 0;
    } else if (isam_has_path_suffix(alias, ".ssam")) {
        alias[strlen(alias) - 5] = 0;
    } else if (isam_has_path_suffix(alias, ".bam")) {
        alias[strlen(alias) - 4] = 0;
    } else if (isam_has_path_suffix(alias, ".sam")) {
        alias[strlen(alias) - 4] = 0;
    } else if (isam_has_path_suffix(alias, ".cram")) {
        alias[strlen(alias) - 5] = 0;
    }
    if (!alias[0]) {
        free(alias);
        alias = xstrdup(type == QUERY_SOURCE_STORE ? "ssam" : "isam");
    }
    return alias;
}

static void check_unique_alias(query_source_t *sources, size_t n,
                               const char *alias) {
    for (size_t i = 0; i < n; ++i)
        if (strcmp(sources[i].alias, alias) == 0)
            die("duplicate query source alias '%s'", alias);
}

static query_source_type_t path_source_type(const char *path) {
    if (isam_is_view_index_path(path))
        return QUERY_SOURCE_VIEW;
    if (isam_is_store_index_path(path) ||
        isam_has_path_suffix(path, ".bam") ||
        isam_has_path_suffix(path, ".cram") ||
        isam_has_path_suffix(path, ".sam"))
        return QUERY_SOURCE_STORE;
    die("cannot infer source type for '%s' (expected .isam, .ssam, BAM, CRAM, or SAM)",
        path);
}

static void add_source(query_source_t **sources,
                       size_t *n_sources,
                       size_t *cap_sources,
                       const char *arg) {
    char *alias = NULL;
    char *path = NULL;
    split_alias_path(arg, &alias, &path);
    query_source_type_t type = path_source_type(path);
    if (!alias)
        alias = default_alias_for_path(path, type);
    check_unique_alias(*sources, *n_sources, alias);

    if (*n_sources == *cap_sources) {
        *cap_sources = *cap_sources ? *cap_sources * 2 : 4;
        *sources = xrealloc(*sources, *cap_sources * sizeof **sources);
    }
    query_source_t *src = &(*sources)[(*n_sources)++];
    memset(src, 0, sizeof *src);
    src->alias = alias;
    src->path = path;
    src->type = type;

    if (type == QUERY_SOURCE_STORE) {
        char *store_path = NULL;
        char *store_index = NULL;
        if (isam_is_store_index_path(path)) {
            lmdb_backend_t *lmdb = lmdb_backend_open(path, 1, 0);
            char *store_type = lmdb_backend_meta_get(lmdb, "isam_type");
            if (!store_type || strcmp(store_type, "store") != 0)
                die("'%s' is not a store .ssam (type=%s)",
                    path, store_type ? store_type : "<missing>");
            free(store_type);
            store_path = lmdb_backend_meta_get(lmdb, "store_path");
            lmdb_backend_close(lmdb);
            if (!store_path)
                die("store index '%s' has no store_path metadata", path);
            store_index = xstrdup(path);
        } else {
            store_path = xstrdup(path);
            store_index = isam_store_index_path_for(path);
        }
        src->store = store_open(store_path, store_index);
        src->hdr = src->store->hdr;
        free(store_path);
        free(store_index);
    } else {
        src->view = view_open_read(path);
        view_link_store(src->view, NULL, NULL);
        src->hdr = src->view->hdr;
    }
}

static int find_source(query_source_t *sources, size_t n_sources,
                       const char *alias) {
    for (size_t i = 0; i < n_sources; ++i)
        if (strcmp(sources[i].alias, alias) == 0)
            return (int)i;
    if (strcmp(alias, "ssam") == 0) {
        int found = -1;
        for (size_t i = 0; i < n_sources; ++i) {
            if (sources[i].type != QUERY_SOURCE_STORE)
                continue;
            if (found >= 0)
                die("alias 'ssam' is ambiguous; use an explicit source alias");
            found = (int)i;
        }
        if (found >= 0)
            return found;
    }
    if (strcmp(alias, "isam") == 0) {
        int found = -1;
        for (size_t i = 0; i < n_sources; ++i) {
            if (sources[i].type != QUERY_SOURCE_VIEW)
                continue;
            if (found >= 0)
                die("alias 'isam' is ambiguous; use an explicit source alias");
            found = (int)i;
        }
        if (found >= 0)
            return found;
    }
    die("unknown query source alias '%s'", alias);
}

static void source_destroy(query_source_t *src) {
    free(src->alias);
    free(src->path);
    if (src->view) {
        view_close(src->view);
    } else if (src->store) {
        store_close(src->store);
    }
}

static char *bam_cigar_string(int n_cigar, const uint32_t *cigar) {
    if (n_cigar <= 0)
        return xstrdup("*");
    size_t cap = (size_t)n_cigar * 12 + 1;
    char *out = xmalloc(cap);
    char *p = out;
    size_t left = cap;
    for (int i = 0; i < n_cigar; ++i) {
        int written = snprintf(p, left, "%u%c",
                               bam_cigar_oplen(cigar[i]),
                               bam_cigar_opchr(cigar[i]));
        if (written < 0 || (size_t)written >= left)
            die("internal CIGAR formatting overflow");
        p += written;
        left -= (size_t)written;
    }
    return out;
}

static void hits_add(query_hits_t *hits,
                     int tid,
                     int64_t pos0,
                     int64_t end0,
                     int flag,
                     int mapq,
                     const char *cigar,
                     int is_unmapped) {
    if (hits->n_alns == hits->cap_alns) {
        hits->cap_alns = hits->cap_alns ? hits->cap_alns * 2 : 2;
        hits->alns = xrealloc(hits->alns, hits->cap_alns * sizeof *hits->alns);
    }
    query_aln_t *aln = &hits->alns[hits->n_alns++];
    memset(aln, 0, sizeof *aln);
    aln->tid = tid;
    aln->pos0 = pos0;
    aln->end0 = end0;
    aln->flag = flag;
    aln->mapq = mapq;
    aln->cigar = xstrdup(cigar ? cigar : ".");
    aln->is_unmapped = is_unmapped;
}

static void hits_destroy(query_hits_t *hits) {
    for (size_t i = 0; i < hits->n_alns; ++i)
        free(hits->alns[i].cigar);
    free(hits->alns);
    memset(hits, 0, sizeof *hits);
}

static void interval_add(query_interval_t **items,
                         size_t *n_items,
                         size_t *cap_items,
                         int tid,
                         int64_t beg,
                         int64_t end) {
    if (tid < 0 || end <= beg)
        return;
    if (*n_items == *cap_items) {
        *cap_items = *cap_items ? *cap_items * 2 : 1024;
        *items = xrealloc(*items, *cap_items * sizeof **items);
    }
    (*items)[*n_items].tid = tid;
    (*items)[*n_items].beg = beg;
    (*items)[*n_items].end = end;
    ++*n_items;
}

static int interval_cmp(const void *a, const void *b) {
    const query_interval_t *x = a;
    const query_interval_t *y = b;
    if (x->tid != y->tid)
        return x->tid < y->tid ? -1 : 1;
    if (x->beg != y->beg)
        return x->beg < y->beg ? -1 : 1;
    if (x->end != y->end)
        return x->end < y->end ? -1 : 1;
    return 0;
}

static void select_region(query_source_t *src,
                          const char *region,
                          qname_set_t *out);

static void collect_source_hits(query_source_t *src,
                                const char *qname,
                                query_hits_t *hits) {
    memset(hits, 0, sizeof *hits);
    if (src->type == QUERY_SOURCE_STORE) {
        store_record_t rec;
        if (!store_lookup_qname(src->store, qname, &rec))
            return;
        hits->has_read = 1;
        hits->read_id = rec.read_id;
        int mapq = -1;
        char *cigar = NULL;
        if (src->needs_payload) {
            bam1_t *payload = bam_init1();
            int rc = store_fetch_payload(src->store, qname, &rec, payload);
            if (rc != 0)
                die("cannot fetch canonical payload for '%s' from '%s' (rc=%d)",
                    qname, src->store->store_path, rc);
            mapq = payload->core.qual;
            cigar = bam_cigar_string(payload->core.n_cigar,
                                     bam_get_cigar(payload));
            bam_destroy1(payload);
        }
        hits_add(hits, rec.tid, rec.pos0, rec.end0, rec.flag, mapq,
                 cigar ? cigar : ".", rec.is_unmapped);
        free(cigar);
        store_record_clear(&rec);
        return;
    }

    int64_t read_id = 0;
    if (!lmdb_view_lookup_qname(src->view->lmdb, qname_hash(qname), qname,
                                &read_id))
        return;
    hits->has_read = 1;
    hits->read_id = read_id;
    lmdb_view_iter_t *it = lmdb_view_iter_open_read_id(src->view->lmdb,
                                                       read_id);
    lmdb_view_record_t rec;
    while (lmdb_view_iter_next(it, &rec)) {
        char *cigar = bam_cigar_string(rec.n_cigar, rec.cigar);
        hits_add(hits, rec.tid, rec.pos0, rec.end0, rec.flag, rec.mapq,
                 cigar, rec.is_unmapped);
        free(cigar);
        lmdb_view_record_clear(&rec);
    }
    lmdb_view_iter_close(it);
}

static void select_interval(query_source_t *src,
                            int tid,
                            int64_t beg,
                            int64_t end,
                            qname_set_t *out) {
    const char *rname = sam_hdr_tid2name(src->hdr, tid);
    if (!rname)
        die("source '%s' has no reference name for tid %d", src->alias, tid);
    char region[1024];
    int n = snprintf(region, sizeof region, "%s:%lld-%lld",
                     rname, (long long)beg + 1, (long long)end);
    if (n < 0 || (size_t)n >= sizeof region)
        die("reference name too long while building query span for '%s'",
            src->alias);
    select_region(src, region, out);
}

static void apply_span(query_source_t *src,
                       qname_set_t *selected,
                       int64_t max_gap) {
    if (max_gap < 0)
        die("--span distance must be non-negative");
    qname_set_sort_unique(selected);

    query_interval_t *intervals = NULL;
    size_t n_intervals = 0, cap_intervals = 0;
    for (size_t i = 0; i < selected->n; ++i) {
        query_hits_t hits;
        collect_source_hits(src, selected->items[i], &hits);
        for (size_t j = 0; j < hits.n_alns; ++j) {
            query_aln_t *aln = &hits.alns[j];
            if (aln->is_unmapped)
                continue;
            interval_add(&intervals, &n_intervals, &cap_intervals,
                         aln->tid, aln->pos0, aln->end0);
        }
        hits_destroy(&hits);
    }

    if (n_intervals == 0) {
        free(intervals);
        return;
    }
    qsort(intervals, n_intervals, sizeof *intervals, interval_cmp);

    qname_set_t expanded = { 0 };
    query_interval_t cur = intervals[0];
    for (size_t i = 1; i < n_intervals; ++i) {
        query_interval_t next = intervals[i];
        if (next.tid == cur.tid &&
            (next.beg <= cur.end || next.beg - cur.end <= max_gap)) {
            if (next.end > cur.end)
                cur.end = next.end;
            continue;
        }
        select_interval(src, cur.tid, cur.beg, cur.end, &expanded);
        cur = next;
    }
    select_interval(src, cur.tid, cur.beg, cur.end, &expanded);

    qname_set_union_into(selected, &expanded);
    qname_set_sort_unique(selected);
    qname_set_destroy(&expanded);
    free(intervals);
}

static int overlaps(int tid, int64_t beg, int64_t end,
                    int rec_tid, int64_t rec_beg, int64_t rec_end) {
    return rec_tid == tid && rec_beg < end && rec_end > beg;
}

static void select_region_store(query_source_t *src,
                                const char *region,
                                qname_set_t *out) {
    hts_idx_t *idx = sam_index_load(src->store->sam, src->store->store_path);
    if (!idx)
        die("cannot load index for store '%s' (samtools index it first)",
            src->store->store_path);
    hts_itr_t *itr = sam_itr_querys(idx, src->store->hdr, region);
    if (!itr)
        die("invalid or unavailable region '%s' for source '%s'",
            region, src->alias);
    bam1_t *record = bam_init1();
    int rc;
    while ((rc = sam_itr_next(src->store->sam, itr, record)) >= 0) {
        const char *qname = bam_get_qname(record);
        store_record_t srec;
        if (store_lookup_qname(src->store, qname, &srec)) {
            if (!srec.is_unmapped &&
                overlaps(srec.tid, srec.pos0, srec.end0,
                         record->core.tid, record->core.pos,
                         bam_endpos(record)))
                qname_set_add(out, qname);
            store_record_clear(&srec);
        }
    }
    bam_destroy1(record);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
}

static void select_region_view(query_source_t *src,
                               const char *region,
                               qname_set_t *out) {
    int64_t has_iv = 0;
    if (!view_meta_get_i64(src->view, "indexed_interval", &has_iv) || !has_iv)
        die("source '%s' has no interval index; run `isamtools sort --interval %s`",
            src->alias, src->view->path);
    int tid = -1;
    hts_pos_t beg = 0, end = HTS_POS_MAX;
    const char *tail = sam_parse_region(src->view->hdr, region, &tid, &beg,
                                        &end, 0);
    if (!tail || tid < 0)
        die("invalid region '%s' for source '%s'", region, src->alias);
    lmdb_view_iter_t *it = lmdb_view_iter_open(src->view->lmdb,
                                               LMDB_VIEW_ITER_REGION,
                                               tid, (int64_t)beg,
                                               (int64_t)end);
    lmdb_view_record_t rec;
    while (lmdb_view_iter_next(it, &rec)) {
        store_record_t srec;
        if (!store_lookup_read_id(src->view->store, rec.store_read_id, &srec))
            die("store_read_id '%lld' missing from linked store for '%s'",
                (long long)rec.store_read_id, src->alias);
        qname_set_add(out, srec.qname);
        store_record_clear(&srec);
        lmdb_view_record_clear(&rec);
    }
    lmdb_view_iter_close(it);
}

static void select_region(query_source_t *src,
                          const char *region,
                          qname_set_t *out) {
    if (src->type == QUERY_SOURCE_STORE)
        select_region_store(src, region, out);
    else
        select_region_view(src, region, out);
}

static int source_subset_lookup(query_source_t *src,
                                const char *name,
                                lmdb_backend_t **lmdb,
                                int64_t *subset_id) {
    if (src->type == QUERY_SOURCE_STORE) {
        if (lmdb_subset_lookup(src->store->lmdb, name, subset_id)) {
            *lmdb = src->store->lmdb;
            return 1;
        }
        return 0;
    }
    if (lmdb_subset_lookup(src->view->lmdb, name, subset_id)) {
        *lmdb = src->view->lmdb;
        return 1;
    }
    if (lmdb_subset_lookup(src->view->store->lmdb, name, subset_id)) {
        *lmdb = src->view->store->lmdb;
        return 1;
    }
    return 0;
}

static void select_subset(query_source_t *src,
                          const char *subset,
                          qname_set_t *out) {
    lmdb_backend_t *subset_lmdb = NULL;
    int64_t subset_id = 0;
    if (!source_subset_lookup(src, subset, &subset_lmdb, &subset_id))
        die("subset '%s' not found in source '%s'", subset, src->alias);

    if (src->type == QUERY_SOURCE_VIEW && subset_lmdb == src->view->lmdb) {
        lmdb_view_iter_t *it = lmdb_view_iter_open_subset(src->view->lmdb,
                                                          subset_id);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            store_record_t srec;
            if (!store_lookup_read_id(src->view->store, rec.store_read_id, &srec))
                die("store_read_id '%lld' missing from linked store for '%s'",
                    (long long)rec.store_read_id, src->alias);
            qname_set_add(out, srec.qname);
            store_record_clear(&srec);
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
        return;
    }

    lmdb_store_iter_t *it = lmdb_store_iter_open(
        src->type == QUERY_SOURCE_STORE ? src->store->lmdb : src->view->store->lmdb);
    lmdb_store_record_t rec;
    while (lmdb_store_iter_next(it, &rec)) {
        if (lmdb_subset_contains_read(subset_lmdb, subset_id, rec.read_id))
            qname_set_add(out, rec.qname);
        lmdb_store_record_clear(&rec);
    }
    lmdb_store_iter_close(it);
}

static void load_qnames_file(const char *path, qname_set_t *out) {
    FILE *fp = strcmp(path, "-") == 0 ? stdin : fopen(path, "r");
    if (!fp)
        die("cannot open qname file '%s': %s", path, strerror(errno));
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        (void)len;
        char *p = line;
        while (*p == ' ' || *p == '\t') ++p;
        if (!*p || *p == '\n' || *p == '#') continue;
        char *end = p;
        while (*end && *end != '\n' && *end != '\r' &&
               *end != '\t' && *end != ' ')
            ++end;
        *end = 0;
        qname_set_add(out, p);
    }
    free(line);
    if (fp != stdin && fclose(fp) != 0)
        die("error closing qname file '%s': %s", path, strerror(errno));
}

static query_field_t parse_field(const char *field) {
    if (strcmp(field, "read_id") == 0) return QUERY_FIELD_READ_ID;
    if (strcmp(field, "n") == 0 || strcmp(field, "n_aln") == 0)
        return QUERY_FIELD_N;
    if (strcmp(field, "rname") == 0) return QUERY_FIELD_RNAME;
    if (strcmp(field, "pos") == 0) return QUERY_FIELD_POS;
    if (strcmp(field, "start0") == 0) return QUERY_FIELD_START0;
    if (strcmp(field, "end0") == 0) return QUERY_FIELD_END0;
    if (strcmp(field, "region") == 0) return QUERY_FIELD_REGION;
    if (strcmp(field, "flag") == 0) return QUERY_FIELD_FLAG;
    if (strcmp(field, "mapq") == 0) return QUERY_FIELD_MAPQ;
    if (strcmp(field, "strand") == 0) return QUERY_FIELD_STRAND;
    if (strcmp(field, "cigar") == 0) return QUERY_FIELD_CIGAR;
    if (strcmp(field, "is_mapped") == 0) return QUERY_FIELD_IS_MAPPED;
    die("unknown query field '%s'", field);
}

static void add_column(query_column_t **cols,
                       size_t *n_cols,
                       size_t *cap_cols,
                       query_source_t *sources,
                       size_t n_sources,
                       const char *spec) {
    if (*n_cols == *cap_cols) {
        *cap_cols = *cap_cols ? *cap_cols * 2 : 16;
        *cols = xrealloc(*cols, *cap_cols * sizeof **cols);
    }
    query_column_t *col = &(*cols)[(*n_cols)++];
    memset(col, 0, sizeof *col);
    col->header = xstrdup(spec);
    if (strcmp(spec, "qname") == 0) {
        col->kind = QUERY_COL_QNAME;
        return;
    }

    const char *dot = strchr(spec, '.');
    if (!dot || dot == spec || !dot[1])
        die("invalid query column '%s' (expected qname or ALIAS.FIELD)", spec);
    char *alias = xstrndup_(spec, (size_t)(dot - spec));
    int source_index = find_source(sources, n_sources, alias);
    free(alias);
    col->source_index = source_index;

    const char *field = dot + 1;
    if (strncmp(field, "in_subset:", 10) == 0) {
        if (!field[10])
            die("invalid query column '%s' (empty subset name)", spec);
        col->kind = QUERY_COL_SOURCE_SUBSET;
        col->subset_name = xstrdup(field + 10);
        return;
    }
    col->kind = QUERY_COL_SOURCE_FIELD;
    col->field = parse_field(field);
    if (col->field == QUERY_FIELD_MAPQ || col->field == QUERY_FIELD_CIGAR)
        sources[source_index].needs_payload = 1;
}

static void parse_columns(const char *arg,
                          query_column_t **cols,
                          size_t *n_cols,
                          size_t *cap_cols,
                          query_source_t *sources,
                          size_t n_sources) {
    char *copy = xstrdup(arg);
    char *save = NULL;
    for (char *tok = strtok_r(copy, ",", &save); tok;
         tok = strtok_r(NULL, ",", &save)) {
        if (!*tok)
            die("empty query column in '%s'", arg);
        add_column(cols, n_cols, cap_cols, sources, n_sources, tok);
    }
    free(copy);
}

static char field_sep(size_t i) {
    return i == 0 ? 0 : ';';
}

static const char *aln_rname(query_source_t *src, const query_aln_t *aln) {
    if (aln->is_unmapped || aln->tid < 0)
        return "*";
    const char *name = sam_hdr_tid2name(src->hdr, aln->tid);
    return name ? name : "<missing>";
}

static void print_source_field(FILE *out,
                               query_source_t *src,
                               const query_hits_t *hits,
                               query_field_t field) {
    if (!hits->has_read) {
        fputc('.', out);
        return;
    }
    if (field == QUERY_FIELD_READ_ID) {
        fprintf(out, "%lld", (long long)hits->read_id);
        return;
    }
    if (field == QUERY_FIELD_N) {
        fprintf(out, "%zu", hits->n_alns);
        return;
    }
    if (hits->n_alns == 0) {
        fputc('.', out);
        return;
    }

    for (size_t i = 0; i < hits->n_alns; ++i) {
        const query_aln_t *aln = &hits->alns[i];
        char sep = field_sep(i);
        if (sep) fputc(sep, out);
        switch (field) {
        case QUERY_FIELD_RNAME:
            fputs(aln_rname(src, aln), out);
            break;
        case QUERY_FIELD_POS:
            if (aln->is_unmapped || aln->pos0 < 0) fputc('.', out);
            else fprintf(out, "%lld", (long long)aln->pos0 + 1);
            break;
        case QUERY_FIELD_START0:
            if (aln->is_unmapped || aln->pos0 < 0) fputc('.', out);
            else fprintf(out, "%lld", (long long)aln->pos0);
            break;
        case QUERY_FIELD_END0:
            if (aln->is_unmapped || aln->end0 < 0) fputc('.', out);
            else fprintf(out, "%lld", (long long)aln->end0);
            break;
        case QUERY_FIELD_REGION:
            if (aln->is_unmapped || aln->pos0 < 0 || aln->end0 <= aln->pos0) {
                fputc('.', out);
            } else {
                fprintf(out, "%s:%lld-%lld",
                        aln_rname(src, aln),
                        (long long)aln->pos0 + 1,
                        (long long)aln->end0);
            }
            break;
        case QUERY_FIELD_FLAG:
            fprintf(out, "%d", aln->flag);
            break;
        case QUERY_FIELD_MAPQ:
            if (aln->mapq < 0) fputc('.', out);
            else fprintf(out, "%d", aln->mapq);
            break;
        case QUERY_FIELD_STRAND:
            fputc((aln->flag & BAM_FREVERSE) ? '-' : '+', out);
            break;
        case QUERY_FIELD_CIGAR:
            fputs(aln->cigar ? aln->cigar : ".", out);
            break;
        case QUERY_FIELD_IS_MAPPED:
            fputc(aln->is_unmapped ? '0' : '1', out);
            break;
        default:
            die("internal query field dispatch error");
        }
    }
}

static void print_subset_value(FILE *out,
                               query_source_t *src,
                               const query_hits_t *hits,
                               const char *subset_name) {
    if (!hits->has_read) {
        fputc('0', out);
        return;
    }
    lmdb_backend_t *lmdb = NULL;
    int64_t subset_id = 0;
    if (!source_subset_lookup(src, subset_name, &lmdb, &subset_id)) {
        fputc('0', out);
        return;
    }
    fputc(lmdb_subset_contains_read(lmdb, subset_id, hits->read_id) ? '1' : '0',
          out);
}

static void column_destroy(query_column_t *col) {
    free(col->header);
    free(col->subset_name);
}

static void selector_set_destroy(selector_set_t *sel) {
    qname_set_destroy(&sel->qnames);
    free(sel->label);
}

static void span_spec_destroy(span_spec_t *span) {
    free(span->arg);
}

int cmd_query(int argc, char **argv) {
    query_source_t *sources = NULL;
    size_t n_sources = 0, cap_sources = 0;
    selector_set_t *selectors = NULL;
    size_t n_selectors = 0, cap_selectors = 0;
    span_spec_t *spans = NULL;
    size_t n_spans = 0, cap_spans = 0;
    query_column_t *columns = NULL;
    size_t n_columns = 0, cap_columns = 0;
    const char *columns_arg = NULL;
    const char *output_path = "-";
    int print_header = 1;
    int join_inner = 0;

    static struct option longopts[] = {
        {"source",    required_argument, 0, 1000},
        {"select",    required_argument, 0, 1001},
        {"region",    required_argument, 0, 1001},
        {"subset",    required_argument, 0, 1002},
        {"qnames",    required_argument, 0, 1003},
        {"columns",   required_argument, 0, 1004},
        {"join",      required_argument, 0, 1005},
        {"no-header", no_argument,       0, 1006},
        {"span",      required_argument, 0, 1007},
        {"output",    required_argument, 0, 'o'},
        {"help",      no_argument,       0, 'h'},
        {0,0,0,0}
    };

    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "o:h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000:
            add_source(&sources, &n_sources, &cap_sources, optarg);
            break;
        case 1001:
        case 1002:
        case 1003:
            if (n_selectors == cap_selectors) {
                cap_selectors = cap_selectors ? cap_selectors * 2 : 8;
                selectors = xrealloc(selectors,
                                     cap_selectors * sizeof *selectors);
            }
            memset(&selectors[n_selectors], 0, sizeof selectors[n_selectors]);
            selectors[n_selectors].label = xstrdup(optarg);
            selectors[n_selectors].qnames.sorted = 0;
            if (c == 1003) {
                load_qnames_file(optarg, &selectors[n_selectors].qnames);
            } else {
                if (n_sources == 0)
                    die("--source must appear before --select/--subset");
                char *alias = NULL, *value = NULL;
                split_alias_value(optarg,
                                  c == 1001 ? "--select" : "--subset",
                                  &alias, &value);
                int src_i = find_source(sources, n_sources, alias);
                if (c == 1001)
                    select_region(&sources[src_i], value,
                                  &selectors[n_selectors].qnames);
                else
                    select_subset(&sources[src_i], value,
                                  &selectors[n_selectors].qnames);
                free(alias);
                free(value);
            }
            ++n_selectors;
            break;
        case 1004:
            columns_arg = optarg;
            break;
        case 1005:
            if (strcmp(optarg, "union") == 0) {
                join_inner = 0;
            } else if (strcmp(optarg, "inner") == 0 ||
                       strcmp(optarg, "intersection") == 0) {
                join_inner = 1;
            } else {
                die("--join must be 'union' or 'inner'");
            }
            break;
        case 1006:
            print_header = 0;
            break;
        case 1007:
            if (n_spans == cap_spans) {
                cap_spans = cap_spans ? cap_spans * 2 : 4;
                spans = xrealloc(spans, cap_spans * sizeof *spans);
            }
            spans[n_spans++].arg = xstrdup(optarg);
            break;
        case 'o':
            output_path = optarg;
            break;
        case 'h':
            usage();
            return 0;
        default:
            usage();
            return 2;
        }
    }
    if (optind < argc) {
        for (; optind < argc; ++optind)
            add_source(&sources, &n_sources, &cap_sources, argv[optind]);
    }
    if (n_sources == 0) {
        usage();
        return 2;
    }
    if (n_selectors == 0)
        die("query requires at least one read selector (--select, --subset, or --qnames)");

    if (!columns_arg) {
        add_column(&columns, &n_columns, &cap_columns,
                   sources, n_sources, "qname");
        for (size_t i = 0; i < n_sources; ++i) {
            char buf[256];
            snprintf(buf, sizeof buf, "%s.region", sources[i].alias);
            add_column(&columns, &n_columns, &cap_columns,
                       sources, n_sources, buf);
        }
    } else {
        parse_columns(columns_arg, &columns, &n_columns, &cap_columns,
                      sources, n_sources);
    }

    qname_set_t selected = { 0 };
    if (join_inner) {
        qname_set_union_into(&selected, &selectors[0].qnames);
        for (size_t i = 1; i < n_selectors; ++i) {
            qname_set_t next = { 0 };
            qname_set_intersect_into(&next, &selected, &selectors[i].qnames);
            qname_set_destroy(&selected);
            selected = next;
        }
    } else {
        for (size_t i = 0; i < n_selectors; ++i)
            qname_set_union_into(&selected, &selectors[i].qnames);
    }
    qname_set_sort_unique(&selected);

    for (size_t i = 0; i < n_spans; ++i) {
        char *alias = NULL, *value = NULL;
        split_alias_value(spans[i].arg, "--span", &alias, &value);
        char *endp = NULL;
        errno = 0;
        long long max_gap = strtoll(value, &endp, 10);
        if (errno || endp == value || *endp || max_gap < 0)
            die("--span distance must be a non-negative integer, got '%s'",
                value);
        int src_i = find_source(sources, n_sources, alias);
        apply_span(&sources[src_i], &selected, (int64_t)max_gap);
        free(alias);
        free(value);
    }

    FILE *out = stdout;
    if (strcmp(output_path, "-") != 0) {
        out = fopen(output_path, "w");
        if (!out)
            die("cannot open output '%s': %s", output_path, strerror(errno));
    }

    if (print_header) {
        for (size_t i = 0; i < n_columns; ++i) {
            if (i) fputc('\t', out);
            fputs(columns[i].header, out);
        }
        fputc('\n', out);
    }

    query_hits_t *hits = xcalloc(n_sources, sizeof *hits);
    for (size_t q = 0; q < selected.n; ++q) {
        for (size_t s = 0; s < n_sources; ++s)
            collect_source_hits(&sources[s], selected.items[q], &hits[s]);

        for (size_t cidx = 0; cidx < n_columns; ++cidx) {
            query_column_t *col = &columns[cidx];
            if (cidx) fputc('\t', out);
            if (col->kind == QUERY_COL_QNAME) {
                fputs(selected.items[q], out);
            } else if (col->kind == QUERY_COL_SOURCE_FIELD) {
                query_source_t *src = &sources[col->source_index];
                print_source_field(out, src, &hits[col->source_index],
                                   col->field);
            } else if (col->kind == QUERY_COL_SOURCE_SUBSET) {
                query_source_t *src = &sources[col->source_index];
                print_subset_value(out, src, &hits[col->source_index],
                                   col->subset_name);
            }
        }
        fputc('\n', out);

        for (size_t s = 0; s < n_sources; ++s)
            hits_destroy(&hits[s]);
    }
    free(hits);

    if (out != stdout && fclose(out) != 0)
        die("error closing output '%s': %s", output_path, strerror(errno));

    qname_set_destroy(&selected);
    for (size_t i = 0; i < n_selectors; ++i)
        selector_set_destroy(&selectors[i]);
    free(selectors);
    for (size_t i = 0; i < n_spans; ++i)
        span_spec_destroy(&spans[i]);
    free(spans);
    for (size_t i = 0; i < n_columns; ++i)
        column_destroy(&columns[i]);
    free(columns);
    for (size_t i = 0; i < n_sources; ++i)
        source_destroy(&sources[i]);
    free(sources);
    return 0;
}
