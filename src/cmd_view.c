#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"
#include "aux_tags.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <strings.h>
#include <sys/types.h>

typedef enum view_output_format {
    VIEW_OUT_SAM = 0,
    VIEW_OUT_BAM = 1,
    VIEW_OUT_CRAM = 2
} view_output_format_t;

static void usage(void) {
    fputs(
        "Usage: isamtools view [options] <view.isam> [region]\n"
        "\n"
        "Emit reconstructed SAM/BAM by combining alignment-view records with\n"
        "the canonical store payload.\n"
        "\n"
        "Output options:\n"
        "  -b, --bam        Output BAM.\n"
        "  -C, --cram       Output CRAM (requires -T).\n"
        "  -1, --fast       Fast BAM compression (and default to BAM).\n"
        "  -u, --uncompressed  Uncompressed BAM output (and default to BAM).\n"
        "  -O, --output-fmt FMT  Output format: sam, bam, or cram.\n"
        "  -h, --with-header  Include @ header lines in SAM output.\n"
        "  -H, --header-only  Print only the header and exit.\n"
        "      --no-header  Suppress SAM header lines (default for SAM).\n"
        "  -c, --count      Count records instead of emitting them.\n"
        "  -o, --output P   Write output to file (default stdout).\n"
        "\n"
        "Input/reconstruction options:\n"
        "      --store P    Override store path.\n"
        "      --store-index P  Override store index path.\n"
        "  -T, --reference P Reference FASTA for CRAM payload decoding.\n"
        "  -@, --threads N  Worker threads for htslib I/O.\n"
        "  -L, --targets-file P  BED file of regions to emit.\n"
        "      --regions-file P  BED file of regions to emit.\n"
        "  -M, --use-index  Accepted for samtools compatibility; indexed paths are used automatically.\n"
        "  -P, --fetch-pairs  Emit all alignments for reads selected by region/BED.\n"
        "      --order ORD  full-scan order: coordinate (default) or name.\n"
        "      --subset S   Only emit reads in named subset (view first, then store).\n"
        "\n"
        "Filtering options:\n"
        "  -f INT           Require all FLAG bits in INT.\n"
        "  -F INT           Exclude records with any FLAG bits in INT.\n"
        "      --rf INT     Require at least one FLAG bit in INT.\n"
        "  -G INT           Exclude records with all FLAG bits in INT.\n"
        "  -q INT           Minimum MAPQ.\n"
        "  -m INT           Minimum query bases consumed by CIGAR.\n"
        "  -r STR           Require read group.\n"
        "  -R P             Read-group allow-list file.\n"
        "  -N P             QNAME allow-list file.\n"
        "  -l STR           Require library (via RG:LB header metadata).\n"
        "  -d TAG[:VALUE]   Require aux tag, optionally with value.\n"
        "  -D TAG:P         Require aux tag whose value is listed in P.\n"
        "\n"
        "General options:\n"
        "  -?, --help       Show this help.\n"
        "  -S               Ignored for samtools compatibility.\n"
        "      --no-PG      Accepted; isamtools does not add @PG lines in view.\n",
        stderr);
}

static void unsupported_option(const char *opt) {
    die("view option '%s' is not implemented by isamtools", opt);
}

static int subset_contains(view_t *v, int subset_from_store,
                           int64_t subset_id, int64_t store_read_id) {
    return lmdb_subset_contains_read(subset_from_store ? v->store->lmdb : v->lmdb,
                                     subset_id, store_read_id);
}

typedef struct str_set {
    char **items;
    size_t n;
    size_t cap;
    int sorted;
} str_set_t;

typedef struct i64_set {
    int64_t *items;
    size_t n;
    size_t cap;
    int sorted;
} i64_set_t;

typedef struct tag_filter {
    char tag[3];
    int has_value;
    char *value;
    str_set_t values;
} tag_filter_t;

typedef struct view_filter {
    int require_all;
    int exclude_any;
    int require_any;
    int exclude_all;
    int min_mapq;
    int min_qlen;
    str_set_t read_groups;
    str_set_t libraries;
    int qname_filter_active;
    i64_set_t qname_ids;
    tag_filter_t *tag_filters;
    size_t n_tag_filters;
    size_t cap_tag_filters;
} view_filter_t;

typedef struct bed_interval {
    int tid;
    int64_t beg;
    int64_t end;
} bed_interval_t;

typedef struct bed_filter {
    bed_interval_t *items;
    size_t n;
    size_t cap;
} bed_filter_t;

typedef struct seen_aln_ids {
    int64_t *items;
    size_t n;
    size_t cap;
} seen_aln_ids_t;

static int bed_interval_cmp(const void *a, const void *b) {
    const bed_interval_t *x = a, *y = b;
    if (x->tid != y->tid) return x->tid < y->tid ? -1 : 1;
    if (x->beg != y->beg) return x->beg < y->beg ? -1 : 1;
    if (x->end != y->end) return x->end < y->end ? -1 : 1;
    return 0;
}

static void bed_filter_add(bed_filter_t *bed_filter, int tid,
                           int64_t beg, int64_t end) {
    if (end <= beg) return;
    if (bed_filter->n == bed_filter->cap) {
        bed_filter->cap = bed_filter->cap ? bed_filter->cap * 2 : 128;
        bed_filter->items = xrealloc(bed_filter->items,
                                     bed_filter->cap * sizeof *bed_filter->items);
    }
    bed_filter->items[bed_filter->n++] = (bed_interval_t){ tid, beg, end };
}

static bed_filter_t bed_filter_load(const char *path, sam_hdr_t *hdr) {
    bed_filter_t bed_filter = { 0 };
    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open BED file '%s': %s", path, strerror(errno));
    char *line = NULL;
    size_t line_cap = 0;
    ssize_t len;
    int64_t line_no = 0;
    while ((len = getline(&line, &line_cap, fp)) >= 0) {
        ++line_no;
        char *p = line;
        while (*p == ' ' || *p == '\t') ++p;
        if (*p == 0 || *p == '\n' || *p == '#') continue;
        if (strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0)
            continue;
        char *chrom = strtok(p, "\t \r\n");
        char *beg_s = strtok(NULL, "\t \r\n");
        char *end_s = strtok(NULL, "\t \r\n");
        if (!chrom || !beg_s || !end_s)
            die("invalid BED line %lld in '%s'", (long long)line_no, path);
        char *beg_end = NULL, *end_end = NULL;
        long long beg = strtoll(beg_s, &beg_end, 10);
        long long end = strtoll(end_s, &end_end, 10);
        if (*beg_end || *end_end || beg < 0 || end < 0 || end < beg)
            die("invalid BED interval on line %lld in '%s'", (long long)line_no, path);
        int tid = sam_hdr_name2tid(hdr, chrom);
        if (tid < 0)
            die("BED reference '%s' on line %lld is absent from the view header",
                chrom, (long long)line_no);
        bed_filter_add(&bed_filter, tid, beg, end);
    }
    free(line);
    fclose(fp);
    qsort(bed_filter.items, bed_filter.n, sizeof *bed_filter.items,
          bed_interval_cmp);
    return bed_filter;
}

static void bed_filter_destroy(bed_filter_t *bed_filter) {
    free(bed_filter->items);
    memset(bed_filter, 0, sizeof *bed_filter);
}

static int overlaps_bed(const bed_filter_t *bed_filter,
                        int tid, int64_t beg, int64_t end) {
    if (!bed_filter) return 1;
    if (bed_filter->n == 0) return 0;
    if (tid < 0 || beg < 0 || end <= beg) return 0;
    for (size_t i = 0; i < bed_filter->n; ++i) {
        const bed_interval_t *iv = &bed_filter->items[i];
        if (iv->tid < tid) continue;
        if (iv->tid > tid) return 0;
        if (iv->beg >= end) return 0;
        if (iv->end > beg) return 1;
    }
    return 0;
}

static int seen_aln_id(seen_aln_ids_t *seen, int64_t aln_id) {
    for (size_t i = 0; i < seen->n; ++i)
        if (seen->items[i] == aln_id) return 1;
    if (seen->n == seen->cap) {
        seen->cap = seen->cap ? seen->cap * 2 : 1024;
        seen->items = xrealloc(seen->items, seen->cap * sizeof *seen->items);
    }
    seen->items[seen->n++] = aln_id;
    return 0;
}

static void seen_aln_ids_destroy(seen_aln_ids_t *seen) {
    free(seen->items);
    memset(seen, 0, sizeof *seen);
}

static int str_ptr_cmp(const void *a, const void *b) {
    const char *const *x = a;
    const char *const *y = b;
    return strcmp(*x, *y);
}

static void str_set_add(str_set_t *set, const char *s) {
    if (set->n == set->cap) {
        set->cap = set->cap ? set->cap * 2 : 32;
        set->items = xrealloc(set->items, set->cap * sizeof *set->items);
    }
    set->items[set->n++] = xstrdup(s);
    set->sorted = 0;
}

static void str_set_sort(str_set_t *set) {
    if (!set->sorted) {
        qsort(set->items, set->n, sizeof *set->items, str_ptr_cmp);
        set->sorted = 1;
    }
}

static int str_set_contains(str_set_t *set, const char *s) {
    if (!set || set->n == 0 || !s) return 0;
    str_set_sort(set);
    char *key = (char *)s;
    return bsearch(&key, set->items, set->n, sizeof *set->items,
                   str_ptr_cmp) != NULL;
}

static void str_set_destroy(str_set_t *set) {
    for (size_t i = 0; i < set->n; ++i)
        free(set->items[i]);
    free(set->items);
    memset(set, 0, sizeof *set);
}

static int i64_cmp(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a;
    int64_t y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static void i64_set_add(i64_set_t *set, int64_t v) {
    if (set->n == set->cap) {
        set->cap = set->cap ? set->cap * 2 : 1024;
        set->items = xrealloc(set->items, set->cap * sizeof *set->items);
    }
    set->items[set->n++] = v;
    set->sorted = 0;
}

static void i64_set_sort(i64_set_t *set) {
    if (!set->sorted) {
        qsort(set->items, set->n, sizeof *set->items, i64_cmp);
        size_t out = 0;
        for (size_t i = 0; i < set->n; ++i)
            if (out == 0 || set->items[i] != set->items[out - 1])
                set->items[out++] = set->items[i];
        set->n = out;
        set->sorted = 1;
    }
}

static int i64_set_contains(i64_set_t *set, int64_t v) {
    if (!set || set->n == 0) return 0;
    i64_set_sort(set);
    return bsearch(&v, set->items, set->n, sizeof *set->items, i64_cmp) != NULL;
}

static void i64_set_destroy(i64_set_t *set) {
    free(set->items);
    memset(set, 0, sizeof *set);
}

static char *first_token(char *line) {
    char *p = line;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
    if (*p == 0 || *p == '#') return NULL;
    char *end = p;
    while (*end && *end != ' ' && *end != '\t' &&
           *end != '\r' && *end != '\n') ++end;
    *end = 0;
    return p;
}

static void str_set_load_file(str_set_t *set, const char *path,
                              const char *what) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open %s file '%s': %s", what, path, strerror(errno));
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        (void)len;
        char *tok = first_token(line);
        if (tok) str_set_add(set, tok);
    }
    free(line);
    if (fclose(fp) != 0)
        die("error closing %s file '%s': %s", what, path, strerror(errno));
    str_set_sort(set);
}

static void view_filter_add_tag(view_filter_t *f, const char *arg) {
    if (strlen(arg) < 2)
        die("invalid tag filter '%s'", arg);
    if (f->n_tag_filters == f->cap_tag_filters) {
        f->cap_tag_filters = f->cap_tag_filters ? f->cap_tag_filters * 2 : 8;
        f->tag_filters = xrealloc(f->tag_filters,
                                  f->cap_tag_filters * sizeof *f->tag_filters);
    }
    tag_filter_t *tf = &f->tag_filters[f->n_tag_filters++];
    memset(tf, 0, sizeof *tf);
    tf->tag[0] = arg[0];
    tf->tag[1] = arg[1];
    tf->tag[2] = 0;
    if (arg[2]) {
        if (arg[2] != ':')
            die("invalid tag filter '%s' (expected TAG or TAG:VALUE)", arg);
        tf->has_value = 1;
        tf->value = xstrdup(arg + 3);
    }
}

static void view_filter_add_tag_file(view_filter_t *f, const char *arg) {
    const char *colon = strchr(arg, ':');
    if (!colon || colon == arg || colon - arg != 2 || !colon[1])
        die("invalid tag-file filter '%s' (expected TAG:FILE)", arg);
    if (f->n_tag_filters == f->cap_tag_filters) {
        f->cap_tag_filters = f->cap_tag_filters ? f->cap_tag_filters * 2 : 8;
        f->tag_filters = xrealloc(f->tag_filters,
                                  f->cap_tag_filters * sizeof *f->tag_filters);
    }
    tag_filter_t *tf = &f->tag_filters[f->n_tag_filters++];
    memset(tf, 0, sizeof *tf);
    tf->tag[0] = arg[0];
    tf->tag[1] = arg[1];
    tf->tag[2] = 0;
    tf->has_value = 1;
    str_set_load_file(&tf->values, colon + 1, "tag value");
}

static void view_filter_destroy(view_filter_t *f) {
    str_set_destroy(&f->read_groups);
    str_set_destroy(&f->libraries);
    i64_set_destroy(&f->qname_ids);
    for (size_t i = 0; i < f->n_tag_filters; ++i) {
        free(f->tag_filters[i].value);
        str_set_destroy(&f->tag_filters[i].values);
    }
    free(f->tag_filters);
    memset(f, 0, sizeof *f);
}

static void load_qname_file(view_t *v, i64_set_t *ids, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("cannot open qname file '%s': %s", path, strerror(errno));
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        (void)len;
        char *qname = first_token(line);
        if (!qname) continue;
        int64_t read_id = 0;
        if (lmdb_view_lookup_qname(v->lmdb, qname_hash(qname), qname, &read_id)) {
            i64_set_add(ids, read_id);
        } else {
            lmdb_store_record_t srec;
            if (lmdb_store_lookup_qname(v->store->lmdb, qname_hash(qname),
                                        qname, &srec)) {
                i64_set_add(ids, srec.read_id);
                lmdb_store_record_clear(&srec);
            }
        }
    }
    free(line);
    if (fclose(fp) != 0)
        die("error closing qname file '%s': %s", path, strerror(errno));
    i64_set_sort(ids);
}

static int parse_flag_arg(const char *s) {
    char *end = NULL;
    long v = strtol(s, &end, 0);
    if (!s[0] || (end && *end) || v < 0 || v > 0xffff)
        die("invalid flag value '%s'", s);
    return (int)v;
}

static int cigar_query_bases(int n_cigar, const uint32_t *cigar) {
    int64_t q = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        if (bam_cigar_type(op) & 1)
            q += bam_cigar_oplen(cigar[i]);
    }
    if (q > INT_MAX) return INT_MAX;
    return (int)q;
}

static int parse_nonnegative_int(const char *s, const char *what) {
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (!s[0] || (end && *end) || v < 0 || v > INT_MAX)
        die("invalid %s '%s'", what, s);
    return (int)v;
}

static void set_output_format(view_output_format_t *fmt,
                              const char *value,
                              int *compression_level) {
    const char *comma = strchr(value, ',');
    size_t n = comma ? (size_t)(comma - value) : strlen(value);
    if (n == 3 && strncasecmp(value, "sam", n) == 0) {
        *fmt = VIEW_OUT_SAM;
    } else if (n == 3 && strncasecmp(value, "bam", n) == 0) {
        *fmt = VIEW_OUT_BAM;
    } else if (n == 4 && strncasecmp(value, "cram", n) == 0) {
        *fmt = VIEW_OUT_CRAM;
    } else {
        die("unknown output format '%.*s' (expected sam, bam, or cram)",
            (int)n, value);
    }

    for (const char *p = comma ? comma + 1 : NULL; p && *p; ) {
        const char *next = strchr(p, ',');
        size_t len = next ? (size_t)(next - p) : strlen(p);
        if (len == 0) {
            /* tolerate repeated commas exactly as an empty option */
        } else if (len == 4 && strncasecmp(p, "fast", len) == 0) {
            *compression_level = 1;
        } else if (len == 12 && strncasecmp(p, "uncompressed", len) == 0) {
            *compression_level = 0;
        } else if (len > 6 && strncasecmp(p, "level=", 6) == 0) {
            char buf[16];
            size_t ncopy = len - 6;
            if (ncopy >= sizeof buf)
                die("invalid output compression option '%.*s'", (int)len, p);
            memcpy(buf, p + 6, ncopy);
            buf[ncopy] = 0;
            int level = parse_nonnegative_int(buf, "compression level");
            if (level > 9)
                die("compression level must be in the range 0..9");
            *compression_level = level;
        } else {
            die("unsupported output format option '%.*s'", (int)len, p);
        }
        p = next ? next + 1 : NULL;
    }
}

static void set_output_format_option(const char *value,
                                     int *compression_level) {
    if (strcasecmp(value, "fast") == 0) {
        *compression_level = 1;
    } else if (strcasecmp(value, "uncompressed") == 0) {
        *compression_level = 0;
    } else if (strncasecmp(value, "level=", 6) == 0) {
        int level = parse_nonnegative_int(value + 6, "compression level");
        if (level > 9)
            die("compression level must be in the range 0..9");
        *compression_level = level;
    } else {
        die("unsupported output format option '%s'", value);
    }
}

static int passes_record_filter(view_filter_t *f,
                                const lmdb_view_record_t *rec) {
    int flag = rec->flag;
    if (f->require_all && (flag & f->require_all) != f->require_all) return 0;
    if (f->exclude_any && (flag & f->exclude_any) != 0) return 0;
    if (f->require_any && (flag & f->require_any) == 0) return 0;
    if (f->exclude_all && (flag & f->exclude_all) == f->exclude_all) return 0;
    if (rec->mapq < f->min_mapq) return 0;
    if (f->min_qlen > 0 &&
        cigar_query_bases(rec->n_cigar, rec->cigar) < f->min_qlen) return 0;
    if (f->qname_filter_active &&
        !i64_set_contains(&f->qname_ids, rec->store_read_id)) return 0;
    return 1;
}

static int view_filter_needs_bam(const view_filter_t *f) {
    return f->read_groups.n > 0 ||
           f->libraries.n > 0 ||
           f->n_tag_filters > 0;
}

static int aux_integer_type(char type) {
    return type == 'c' || type == 'C' || type == 's' || type == 'S' ||
           type == 'i' || type == 'I';
}

static int aux_value_equals(const uint8_t *field, const char *wanted) {
    char type = bam_aux_type(field);
    if (type == 'Z' || type == 'H') {
        const char *s = bam_aux2Z(field);
        return s && strcmp(s, wanted) == 0;
    }
    if (type == 'A') {
        return wanted[0] && !wanted[1] && bam_aux2A(field) == wanted[0];
    }
    if (aux_integer_type(type)) {
        char *end = NULL;
        long long want = strtoll(wanted, &end, 10);
        if (!wanted[0] || (end && *end)) return 0;
        return (long long)bam_aux2i(field) == want;
    }
    if (type == 'f' || type == 'd') {
        char *end = NULL;
        double want = strtod(wanted, &end);
        if (!wanted[0] || (end && *end)) return 0;
        double got = bam_aux2f(field);
        double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
        return fabs(got - want) <= 1e-6 * scale;
    }
    return 0;
}

static int tag_filter_matches(const tag_filter_t *tf, const bam1_t *b) {
    uint8_t *field = bam_aux_get((bam1_t *)b, tf->tag);
    if (!field) return 0;
    if (!tf->has_value) return 1;
    if (tf->value) return aux_value_equals(field, tf->value);
    for (size_t i = 0; i < tf->values.n; ++i)
        if (aux_value_equals(field, tf->values.items[i])) return 1;
    return 0;
}

static int rg_library_allowed(sam_hdr_t *hdr, const char *rg,
                              str_set_t *libraries) {
    if (libraries->n == 0) return 1;
    if (!rg) return 0;
    kstring_t ks = { 0, 0, NULL };
    int rc = sam_hdr_find_tag_id(hdr, "RG", "ID", rg, "LB", &ks);
    int ok = rc == 0 && str_set_contains(libraries, ks.s);
    free(ks.s);
    return ok;
}

static int passes_bam_filter(view_filter_t *f, sam_hdr_t *hdr,
                             const bam1_t *b) {
    if (f->read_groups.n > 0 || f->libraries.n > 0) {
        uint8_t *rg_field = bam_aux_get((bam1_t *)b, "RG");
        const char *rg = rg_field ? bam_aux2Z(rg_field) : NULL;
        if (f->read_groups.n > 0 && !str_set_contains(&f->read_groups, rg))
            return 0;
        if (!rg_library_allowed(hdr, rg, &f->libraries))
            return 0;
    }
    for (size_t i = 0; i < f->n_tag_filters; ++i)
        if (!tag_filter_matches(&f->tag_filters[i], b)) return 0;
    return 1;
}

static void rebuild_lmdb_record(view_t *v, const lmdb_view_record_t *rec,
                                payload_cache_t *cache,
                                bam1_t *canon_buf, bam1_t *out_buf,
                                const isam_tag_set_t *shared_tags,
                                const isam_tag_set_t *dropped_tags) {
    store_record_t srec;
    if (!store_lookup_read_id(v->store, rec->store_read_id, &srec))
        die("store_read_id '%lld' missing from store '%s' (run isamtools rehint)",
            (long long)rec->store_read_id, v->store->store_path);
    int rc = payload_cache_fetch(cache, v->store, &srec, canon_buf);
    if (rc != 0) die("payload fetch failed for '%s' (rc=%d)", srec.qname, rc);

    isam_rebuild_bam(out_buf, canon_buf,
                     rec->flag, rec->tid, rec->pos0, rec->mapq,
                     rec->n_cigar, rec->cigar,
                     rec->mate_tid, rec->mate_pos0, rec->tlen,
                     rec->aux, rec->aux_len);
    lmdb_read_aux_apply(v->lmdb, rec->store_read_id, out_buf, 0);
    lmdb_read_tags_apply(v->lmdb, rec->store_read_id, out_buf, 0);
    isam_aux_tags_append_missing_from(canon_buf, out_buf, shared_tags, dropped_tags);
    lmdb_read_tags_apply(v->store->lmdb, rec->store_read_id, out_buf, 0);
    isam_aux_tags_remove(out_buf, dropped_tags);
    store_record_clear(&srec);
}

static int emit_lmdb_record(samFile *out_sam, sam_hdr_t *hdr,
                            view_t *v, const lmdb_view_record_t *rec,
                            payload_cache_t *cache,
                            bam1_t *canon_buf, bam1_t *out_buf,
                            view_filter_t *filter,
                            const bed_filter_t *bed_filter,
                            const isam_tag_set_t *shared_tags,
                            const isam_tag_set_t *dropped_tags,
                            int count_only) {
    if (!passes_record_filter(filter, rec)) return 0;
    if (!overlaps_bed(bed_filter, rec->tid, rec->pos0, rec->end0)) return 0;
    int need_bam = !count_only || view_filter_needs_bam(filter);
    if (need_bam) {
        rebuild_lmdb_record(v, rec, cache, canon_buf, out_buf,
                            shared_tags, dropped_tags);
        if (!passes_bam_filter(filter, hdr, out_buf)) return 0;
    }
    if (count_only) return 1;

    if (sam_write1(out_sam, hdr, out_buf) < 0)
        die("sam_write1 failed");
    return 1;
}

int cmd_view(int argc, char **argv) {
    const char *store_override = NULL;
    const char *store_idx_override = NULL;
    const char *reference_path = NULL;
    const char *output_path = "-";
    const char *subset = NULL;
    const char *bed_path = NULL;
    const char *qname_path = NULL;
    view_output_format_t out_fmt = VIEW_OUT_SAM;
    int compression_level = -1;
    int print_header = 0;
    int header_only  = 0;
    int count_only = 0;
    int fetch_pairs = 0;
    int n_threads = 0;
    const char *order = "coordinate";
    view_filter_t filter = { 0 };

    static struct option longopts[] = {
        {"bam",         no_argument, 0, 'b'},
        {"cram",        no_argument, 0, 'C'},
        {"fast",        no_argument, 0, '1'},
        {"uncompressed", no_argument, 0, 'u'},
        {"header",      no_argument, 0, 'h'},
        {"with-header", no_argument, 0, 'h'},
        {"header-only", no_argument, 0, 'H'},
        {"no-header",   no_argument, 0, 1000},
        {"store",       required_argument, 0, 1001},
        {"store-index", required_argument, 0, 1002},
        {"order",       required_argument, 0, 1003},
        {"output",      required_argument, 0, 'o'},
        {"output-fmt",  required_argument, 0, 'O'},
        {"output-fmt-option", required_argument, 0, 1007},
        {"help",        no_argument, 0, 1004},
        {"rf",          required_argument, 0, 1005},
        {"incl-flags",  required_argument, 0, 1005},
        {"include-flags", required_argument, 0, 1005},
        {"subset",      required_argument, 0, 1006},
        {"targets-file", required_argument, 0, 'L'},
        {"target-file", required_argument, 0, 'L'},
        {"regions-file", required_argument, 0, 'L'},
        {"count",       no_argument,       0, 'c'},
        {"reference",   required_argument, 0, 'T'},
        {"threads",     required_argument, 0, '@'},
        {"use-index",   no_argument,       0, 'M'},
        {"no-PG",       no_argument,       0, 1008},
        {"unoutput",    required_argument, 0, 'U'},
        {"output-unselected", required_argument, 0, 'U'},
        {"unmap",       no_argument,       0, 'p'},
        {"fetch-pairs", no_argument,       0, 'P'},
        {"fai-reference", required_argument, 0, 't'},
        {"customized-index", no_argument, 0, 'X'},
        {"read-group",  required_argument, 0, 'r'},
        {"read-group-file", required_argument, 0, 'R'},
        {"qname-file",  required_argument, 0, 'N'},
        {"tag",         required_argument, 0, 'd'},
        {"tag-file",    required_argument, 0, 'D'},
        {"library",     required_argument, 0, 'l'},
        {"min-qlen",    required_argument, 0, 'm'},
        {"expr",        required_argument, 0, 'e'},
        {"subsample",   required_argument, 0, 1009},
        {"subsample-seed", required_argument, 0, 1010},
        {"add-flags",   required_argument, 0, 1011},
        {"remove-flags", required_argument, 0, 1012},
        {"remove-tag",  required_argument, 0, 'x'},
        {"keep-tag",    required_argument, 0, 1013},
        {"remove-B",    no_argument,       0, 'B'},
        {"sanitize",    required_argument, 0, 'z'},
        {"write-index", no_argument,       0, 1014},
        {"verbosity",   required_argument, 0, 1015},
        {"input-fmt-option", required_argument, 0, 1016},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv,
                            "1bBCcd:D:e:f:F:G:hHL:l:Mm:N:O:o:PpR:r:Ss:t:T:uU:Xx:z:@:?",
                            longopts, NULL)) != -1) {
        switch (c) {
        case '1': out_fmt = VIEW_OUT_BAM; compression_level = 1; break;
        case 'b': out_fmt = VIEW_OUT_BAM; break;
        case 'B': unsupported_option("-B/--remove-B"); break;
        case 'C': out_fmt = VIEW_OUT_CRAM; break;
        case 'c': count_only = 1; print_header = 0; break;
        case 'd': view_filter_add_tag(&filter, optarg); break;
        case 'D': view_filter_add_tag_file(&filter, optarg); break;
        case 'e': unsupported_option("-e/--expr"); break;
        case 'f': filter.require_all = parse_flag_arg(optarg); break;
        case 'F': filter.exclude_any = parse_flag_arg(optarg); break;
        case 'G': filter.exclude_all = parse_flag_arg(optarg); break;
        case 'h': print_header = 1; break;
        case 'H': header_only = 1; print_header = 1; break;
        case 'L': bed_path = optarg; break;
        case 'l': str_set_add(&filter.libraries, optarg); break;
        case 'M': break;
        case 'm': filter.min_qlen = parse_nonnegative_int(optarg, "minimum query length"); break;
        case 'N': qname_path = optarg; break;
        case 'O':
            set_output_format(&out_fmt, optarg, &compression_level);
            break;
        case 'P': fetch_pairs = 1; break;
        case 'p': unsupported_option("-p/--unmap"); break;
        case 'q': filter.min_mapq = atoi(optarg); break;
        case 'R': str_set_load_file(&filter.read_groups, optarg, "read-group"); break;
        case 'r': str_set_add(&filter.read_groups, optarg); break;
        case 'S': break;
        case 's': unsupported_option("-s/--subsample"); break;
        case 't': break;
        case 'T': reference_path = optarg; break;
        case 'u': out_fmt = VIEW_OUT_BAM; compression_level = 0; break;
        case 'U': unsupported_option("-U/--output-unselected"); break;
        case 'X': unsupported_option("-X/--customized-index"); break;
        case 'x': unsupported_option("-x/--remove-tag"); break;
        case 'z': unsupported_option("-z/--sanitize"); break;
        case '@': n_threads = atoi(optarg); break;
        case '?':
            if (optind > 0 && strcmp(argv[optind - 1], "-?") == 0) {
                usage();
                return 0;
            }
            usage();
            return 2;
        case 1000: print_header = 0; break;
        case 1001: store_override = optarg; break;
        case 1002: store_idx_override = optarg; break;
        case 1003: order = optarg; break;
        case 'o': output_path = optarg; break;
        case 1004: usage(); return 0;
        case 1005: filter.require_any = parse_flag_arg(optarg); break;
        case 1006: subset = optarg; break;
        case 1007: set_output_format_option(optarg, &compression_level); break;
        case 1008: break;
        case 1009: unsupported_option("--subsample"); break;
        case 1010: unsupported_option("--subsample-seed"); break;
        case 1011: unsupported_option("--add-flags"); break;
        case 1012: unsupported_option("--remove-flags"); break;
        case 1013: unsupported_option("--keep-tag"); break;
        case 1014: unsupported_option("--write-index"); break;
        case 1015: break;
        case 1016: unsupported_option("--input-fmt-option"); break;
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *isam_path = argv[optind++];
    int n_regions = argc - optind;
    char **regions = argv + optind;
    const char *region = n_regions == 1 ? regions[0] : NULL;
    if (out_fmt == VIEW_OUT_CRAM && !reference_path)
        die("CRAM output requires -T/--reference");
    if (out_fmt == VIEW_OUT_SAM && compression_level >= 0)
        die("compression options require BAM or CRAM output");
    if (fetch_pairs && n_regions == 0 && !bed_path)
        die("-P/--fetch-pairs requires a region or -L/--targets-file");

    view_t *v = view_open_read(isam_path);
    view_link_store(v, store_override, store_idx_override);
    if (qname_path) {
        filter.qname_filter_active = 1;
        load_qname_file(v, &filter.qname_ids, qname_path);
    }
    bed_filter_t bed_filter = { 0 };
    int use_bed_filter = 0;
    if (bed_path) {
        bed_filter = bed_filter_load(bed_path, v->hdr);
        use_bed_filter = 1;
    }
    if (n_regions > 1) {
        bed_filter_t cli_regions = { 0 };
        for (int i = 0; i < n_regions; ++i) {
            if (strcmp(regions[i], "*") == 0)
                die("'*' unmapped region cannot be combined with other regions");
            int tid = -1;
            hts_pos_t hbeg = 0, hend = HTS_POS_MAX;
            const char *tail = sam_parse_region(v->hdr, regions[i], &tid,
                                                &hbeg, &hend, 0);
            if (!tail || tid < 0) die("invalid region: '%s'", regions[i]);
            bed_filter_add(&cli_regions, tid, (int64_t)hbeg, (int64_t)hend);
        }
        if (use_bed_filter) {
            bed_filter_t intersected = { 0 };
            for (size_t i = 0; i < bed_filter.n; ++i) {
                for (size_t j = 0; j < cli_regions.n; ++j) {
                    if (bed_filter.items[i].tid != cli_regions.items[j].tid)
                        continue;
                    int64_t beg = bed_filter.items[i].beg > cli_regions.items[j].beg
                        ? bed_filter.items[i].beg : cli_regions.items[j].beg;
                    int64_t end = bed_filter.items[i].end < cli_regions.items[j].end
                        ? bed_filter.items[i].end : cli_regions.items[j].end;
                    bed_filter_add(&intersected, bed_filter.items[i].tid, beg, end);
                }
            }
            bed_filter_destroy(&bed_filter);
            bed_filter_destroy(&cli_regions);
            bed_filter = intersected;
        } else {
            bed_filter = cli_regions;
            use_bed_filter = 1;
        }
        qsort(bed_filter.items, bed_filter.n, sizeof *bed_filter.items,
              bed_interval_cmp);
    }
    if (reference_path &&
        hts_set_fai_filename(v->store->sam, reference_path) != 0)
        die("cannot set reference FASTA '%s'", reference_path);
    if (n_threads > 1 &&
        hts_set_threads(v->store->sam, n_threads) != 0)
        die("cannot set store threads to %d", n_threads);
    int have_subset = 0;
    int subset_from_store = 0;
    int64_t subset_id = 0;
    if (subset) {
        if (lmdb_subset_lookup(v->lmdb, subset, &subset_id)) {
            have_subset = 1;
        } else if (lmdb_subset_lookup(v->store->lmdb, subset, &subset_id)) {
            have_subset = 1;
            subset_from_store = 1;
        } else {
            die("subset '%s' is empty or absent in both view and store", subset);
        }
    }

    /* Open output. */
    char mode_str[8];
    char *mp = mode_str;
    *mp++ = 'w';
    if (out_fmt == VIEW_OUT_BAM) *mp++ = 'b';
    else if (out_fmt == VIEW_OUT_CRAM) *mp++ = 'c';
    if (compression_level >= 0) *mp++ = (char)('0' + compression_level);
    *mp = 0;
    samFile *out = NULL;
    if (!count_only) {
        out = sam_open(output_path, mode_str);
        if (!out) die("cannot open output '%s'", output_path);
    }
    if (out && reference_path && out_fmt == VIEW_OUT_CRAM &&
        hts_set_fai_filename(out, reference_path) != 0)
        die("cannot set output reference FASTA '%s'", reference_path);
    if (out && n_threads > 1 && hts_set_threads(out, n_threads) != 0)
        die("cannot set output threads to %d", n_threads);

    if (out && (print_header || out_fmt != VIEW_OUT_SAM)) {
        if (sam_hdr_write(out, v->hdr) < 0) die("sam_hdr_write failed");
    }
    if (header_only) {
        bed_filter_destroy(&bed_filter);
        view_filter_destroy(&filter);
        if (out) sam_close(out);
        view_close(v);
        return 0;
    }

    bam1_t *canon = bam_init1();
    bam1_t *outb  = bam_init1();
    payload_cache_t *cache = payload_cache_create(64 * 1024 * 1024);
    isam_tag_set_t shared_tags;
    isam_tag_set_t dropped_tags;
    isam_tag_set_init(&shared_tags);
    isam_tag_set_init(&dropped_tags);
    char *shared_csv = view_meta_get(v, "shared_tags_snapshot");
    char *dropped_csv = view_meta_get(v, "dropped_tags");
    if (isam_tag_set_parse_csv(&shared_tags, shared_csv) != 0)
        die("invalid shared tag metadata in '%s'", isam_path);
    if (isam_tag_set_parse_csv(&dropped_tags, dropped_csv) != 0)
        die("invalid dropped tag metadata in '%s'", isam_path);
    free(shared_csv);
    free(dropped_csv);
    int64_t n_emit = 0;
    i64_set_t pair_reads = { 0 };

    lmdb_view_iter_kind_t iter_kind = LMDB_VIEW_ITER_IMPORT;
    int tid = -1;
    int64_t beg = 0, end = 0;
    int unmapped_region = 0;
    if (region) {
        if (strcmp(region, "*") == 0) {
            int64_t has_coord = 0;
            view_meta_get_i64(v, "indexed_coord", &has_coord);
            if (!has_coord)
                die("alignment view has no LMDB coordinate index; "
                    "run `isamtools sort %s`", isam_path);
            if (use_bed_filter)
                die("cannot combine '*' unmapped region with region-file filters");
            unmapped_region = 1;
            iter_kind = LMDB_VIEW_ITER_UNMAPPED;
        } else {
            hts_pos_t hbeg = 0, hend = HTS_POS_MAX;
            const char *tail = sam_parse_region(v->hdr, region, &tid,
                                                &hbeg, &hend, 0);
            if (!tail || tid < 0) die("invalid region: '%s'", region);
            int64_t has_iv = 0;
            view_meta_get_i64(v, "indexed_interval", &has_iv);
            if (!has_iv)
                die("alignment view has no LMDB interval index; "
                    "run `isamtools sort --interval %s`", isam_path);
            beg = (int64_t)hbeg;
            end = (int64_t)hend;
            iter_kind = LMDB_VIEW_ITER_REGION;
        }
    }

    if (have_subset && !subset_from_store && (use_bed_filter || region)) {
        lmdb_view_iter_t *it = lmdb_view_iter_open_subset(v->lmdb, subset_id);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            if (unmapped_region && !rec.is_unmapped) {
                lmdb_view_record_clear(&rec);
                continue;
            }
            if (region && !unmapped_region &&
                !(rec.tid == tid && rec.pos0 < end && rec.end0 > beg)) {
                lmdb_view_record_clear(&rec);
                continue;
            }
            int selected = emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                            canon, outb, &filter,
                                            use_bed_filter ? &bed_filter : NULL,
                                            &shared_tags, &dropped_tags,
                                            fetch_pairs ? 1 : count_only);
            if (fetch_pairs) {
                if (selected) i64_set_add(&pair_reads, rec.store_read_id);
            } else {
                n_emit += selected;
            }
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
    } else if (use_bed_filter) {
        int64_t has_iv = 0;
        view_meta_get_i64(v, "indexed_interval", &has_iv);
        if (!has_iv)
            die("alignment view has no LMDB interval index; "
                "run `isamtools sort --interval %s`", isam_path);
        seen_aln_ids_t seen = { 0 };
        for (size_t i = 0; i < bed_filter.n; ++i) {
            int iv_tid = bed_filter.items[i].tid;
            int64_t iv_beg = bed_filter.items[i].beg;
            int64_t iv_end = bed_filter.items[i].end;
            if (region) {
                if (iv_tid != tid || iv_end <= beg || iv_beg >= end)
                    continue;
                if (iv_beg < beg) iv_beg = beg;
                if (iv_end > end) iv_end = end;
            }
            lmdb_view_iter_t *it = lmdb_view_iter_open(v->lmdb,
                                                       LMDB_VIEW_ITER_REGION,
                                                       iv_tid, iv_beg, iv_end);
            lmdb_view_record_t rec;
            while (lmdb_view_iter_next(it, &rec)) {
                if (have_subset && subset_from_store &&
                    !subset_contains(v, subset_from_store, subset_id,
                                     rec.store_read_id)) {
                    lmdb_view_record_clear(&rec);
                    continue;
                }
                if (!seen_aln_id(&seen, rec.aln_id)) {
                    int selected = emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                                    canon, outb, &filter,
                                                    &bed_filter,
                                                    &shared_tags, &dropped_tags,
                                                    fetch_pairs ? 1 : count_only);
                    if (fetch_pairs) {
                        if (selected) i64_set_add(&pair_reads, rec.store_read_id);
                    } else {
                        n_emit += selected;
                    }
                }
                lmdb_view_record_clear(&rec);
            }
            lmdb_view_iter_close(it);
        }
        seen_aln_ids_destroy(&seen);
    } else {
        if (!region && have_subset && !subset_from_store) {
            iter_kind = LMDB_VIEW_ITER_SUBSET;
        } else if (!region && strcasecmp(order, "name") == 0) {
            iter_kind = LMDB_VIEW_ITER_NAME;
        } else if (!region) {
            int64_t has_coord = 0;
            view_meta_get_i64(v, "indexed_coord", &has_coord);
            iter_kind = has_coord ? LMDB_VIEW_ITER_COORD : LMDB_VIEW_ITER_IMPORT;
        }

        lmdb_view_iter_t *it = have_subset && !region
            ? (subset_from_store
                ? lmdb_view_iter_open(v->lmdb, iter_kind, tid, beg, end)
                : lmdb_view_iter_open_subset(v->lmdb, subset_id))
            : lmdb_view_iter_open(v->lmdb, iter_kind, tid, beg, end);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            if (have_subset && subset_from_store &&
                !subset_contains(v, subset_from_store, subset_id,
                                 rec.store_read_id)) {
                lmdb_view_record_clear(&rec);
                continue;
            }
            int selected = emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                            canon, outb, &filter,
                                            use_bed_filter ? &bed_filter : NULL,
                                            &shared_tags, &dropped_tags,
                                            fetch_pairs ? 1 : count_only);
            if (fetch_pairs) {
                if (selected) i64_set_add(&pair_reads, rec.store_read_id);
            } else {
                n_emit += selected;
            }
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
    }

    if (fetch_pairs) {
        i64_set_sort(&pair_reads);
        seen_aln_ids_t emitted = { 0 };
        for (size_t i = 0; i < pair_reads.n; ++i) {
            lmdb_view_iter_t *it = lmdb_view_iter_open_read_id(v->lmdb,
                                                               pair_reads.items[i]);
            lmdb_view_record_t rec;
            while (lmdb_view_iter_next(it, &rec)) {
                if (!seen_aln_id(&emitted, rec.aln_id))
                    n_emit += emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                               canon, outb, &filter, NULL,
                                               &shared_tags, &dropped_tags,
                                               count_only);
                lmdb_view_record_clear(&rec);
            }
            lmdb_view_iter_close(it);
        }
        seen_aln_ids_destroy(&emitted);
    }

    if (count_only) {
        FILE *count_out = stdout;
        if (strcmp(output_path, "-") != 0) {
            count_out = fopen(output_path, "w");
            if (!count_out)
                die("cannot open output '%s': %s", output_path, strerror(errno));
        }
        fprintf(count_out, "%lld\n", (long long)n_emit);
        if (count_out != stdout && fclose(count_out) != 0)
            die("error closing output '%s': %s", output_path, strerror(errno));
    }
    i64_set_destroy(&pair_reads);
    bed_filter_destroy(&bed_filter);
    view_filter_destroy(&filter);
    payload_cache_destroy(cache);
    isam_tag_set_destroy(&shared_tags);
    isam_tag_set_destroy(&dropped_tags);
    bam_destroy1(canon);
    bam_destroy1(outb);
    if (out) sam_close(out);
    view_close(v);
    return 0;
}
