#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"
#include "reference.h"
#include "aux_tags.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <getopt.h>
#include <sys/stat.h>
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
        "      --write-index  Build OUTPUT.bai while writing a coordinate-ordered BAM.\n"
        "\n"
        "Input/reconstruction options:\n"
        "      --store P    Override store path.\n"
        "      --store-index P  Override store index path.\n"
        "  -T, --reference P Reference FASTA for CRAM payload decoding.\n"
        "      --reference-path DIR  Verified, raw M5 reference files (alternative to -T).\n"
        "      --prepare-reference DIR  Prepare a new raw-M5 directory from -T, then use it.\n"
        "  -@, --threads N  Worker threads for htslib I/O.\n"
        "      --decode-threads N  Independent CRAM container decoders (1..64).\n"
        "                   Coordinate output only, whole view or one region.\n"
        "      --decode-lookahead N  Coordinate metadata window (1..65536; default 512).\n"
        "      --cache-selected  Plan requested payloads; release them after final use.\n"
        "                   Requires --decode-threads and an extra metadata pass.\n"
        "  -L, --targets-file P  BED file of regions to emit.\n"
        "      --regions-file P  BED file of regions to emit.\n"
        "  -M, --use-index  Accepted for samtools compatibility; indexed paths are used automatically.\n"
        "  -P, --fetch-pairs  Emit all alignments for reads selected by region/BED.\n"
        "      --order ORD  full-scan order: coordinate (default), name, or store.\n"
        "                   store streams the canonical input once; sort output before indexing.\n"
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
        "      --subsample F  Keep fraction F of read names.\n"
        "      --subsample-seed INT  Seed for --subsample.\n"
        "  -s INT.FRAC      Same as --subsample 0.FRAC --subsample-seed INT.\n"
        "\n"
        "Processing options:\n"
        "  -x TAGS          Comma-separated aux tags to strip; -x ^TAGS keeps only TAGS.\n"
        "      --keep-tag TAGS  Preserve only comma-separated aux tags.\n"
        "  -p, --unmap      Emit unselected records with UNMAP flag set.\n"
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

typedef struct view_output_policy {
    int unmap_unselected;
    isam_tag_set_t remove_tags;
    isam_tag_set_t keep_tags;
    int keep_tags_active;
} view_output_policy_t;

typedef struct view_filter {
    int require_all;
    int exclude_any;
    int require_any;
    int exclude_all;
    int min_mapq;
    int min_qlen;
    int subsample_active;
    double subsample_fraction;
    uint64_t subsample_seed;
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

static void parse_tag_list_into(isam_tag_set_t *tags, const char *arg) {
    if (isam_tag_set_parse_csv(tags, arg) != 0)
        die("invalid tag list '%s'", arg);
}

static void remove_tags_not_in(bam1_t *record, const isam_tag_set_t *keep_tags) {
    int changed;
    do {
        changed = 0;
        for (uint8_t *field = bam_aux_first(record); field;
             field = bam_aux_next(record, field)) {
            char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
            if (isam_tag_set_contains(keep_tags, tag))
                continue;
            if (bam_aux_del(record, field) != 0)
                die("failed to remove aux tag '%s'", tag);
            changed = 1;
            break;
        }
    } while (changed);
}

static void apply_output_policy(bam1_t *record,
                                const view_output_policy_t *policy,
                                int selected) {
    if (!selected && policy->unmap_unselected)
        record->core.flag |= BAM_FUNMAP;
    if (policy->keep_tags_active)
        remove_tags_not_in(record, &policy->keep_tags);
    isam_aux_tags_remove(record, &policy->remove_tags);
}

static void view_output_policy_destroy(view_output_policy_t *policy) {
    isam_tag_set_destroy(&policy->remove_tags);
    isam_tag_set_destroy(&policy->keep_tags);
    memset(policy, 0, sizeof *policy);
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

static uint64_t mix_u64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static int subsample_keep(const view_filter_t *f, int64_t read_id) {
    if (!f->subsample_active)
        return 1;
    if (f->subsample_fraction <= 0.0)
        return 0;
    if (f->subsample_fraction >= 1.0)
        return 1;
    uint64_t h = mix_u64((uint64_t)read_id ^ mix_u64(f->subsample_seed));
    double unit = (double)(h >> 11) * (1.0 / 9007199254740992.0);
    return unit < f->subsample_fraction;
}

static double parse_fraction(const char *s, const char *what) {
    char *end = NULL;
    double v = strtod(s, &end);
    if (!s[0] || (end && *end) || v < 0.0 || v > 1.0)
        die("invalid %s '%s' (expected 0..1)", what, s);
    return v;
}

static int parse_nonnegative_int(const char *s, const char *what);

static void parse_subsample_s(view_filter_t *f, const char *arg) {
    const char *dot = strchr(arg, '.');
    if (!dot || dot == arg || !dot[1])
        die("invalid -s argument '%s' (expected INT.FRAC)", arg);
    char seed_buf[32];
    size_t seed_len = (size_t)(dot - arg);
    if (seed_len >= sizeof seed_buf)
        die("invalid -s seed in '%s'", arg);
    memcpy(seed_buf, arg, seed_len);
    seed_buf[seed_len] = 0;
    f->subsample_seed = (uint64_t)parse_nonnegative_int(seed_buf, "subsample seed");

    char frac_buf[64];
    size_t frac_len = strlen(dot + 1);
    if (frac_len + 3 > sizeof frac_buf)
        die("invalid -s fraction in '%s'", arg);
    frac_buf[0] = '0';
    frac_buf[1] = '.';
    memcpy(frac_buf + 2, dot + 1, frac_len + 1);
    f->subsample_fraction = parse_fraction(frac_buf, "subsample fraction");
    f->subsample_active = 1;
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
    if (!subsample_keep(f, rec->store_read_id)) return 0;
    return 1;
}

static int view_filter_needs_bam(const view_filter_t *f) {
    return f->read_groups.n > 0 ||
           f->libraries.n > 0 ||
           f->n_tag_filters > 0;
}

static int view_filter_is_noop(const view_filter_t *f) {
    return f->require_all == 0 &&
           f->exclude_any == 0 &&
           f->require_any == 0 &&
           f->exclude_all == 0 &&
           f->min_mapq == 0 &&
           f->min_qlen == 0 &&
           !f->subsample_active &&
           f->read_groups.n == 0 &&
           f->libraries.n == 0 &&
           !f->qname_filter_active &&
           f->n_tag_filters == 0;
}

static void write_count(const char *output_path, int64_t count) {
    FILE *count_out = stdout;
    if (strcmp(output_path, "-") != 0) {
        count_out = fopen(output_path, "w");
        if (!count_out)
            die("cannot open output '%s': %s", output_path, strerror(errno));
    }
    fprintf(count_out, "%lld\n", (long long)count);
    if (count_out != stdout && fclose(count_out) != 0)
        die("error closing output '%s': %s", output_path, strerror(errno));
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

static void rebuild_from_payload(view_t *v, const lmdb_view_record_t *rec,
                                const bam1_t *canon_buf, bam1_t *out_buf,
                                const isam_tag_set_t *shared_tags,
                                const isam_tag_set_t *dropped_tags) {
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
}

static int emit_lmdb_record(samFile *out_sam, sam_hdr_t *hdr,
                            view_t *v, const lmdb_view_record_t *rec,
                            payload_cache_t *cache,
                            bam1_t *canon_buf, bam1_t *out_buf,
                            view_filter_t *filter,
                            const bed_filter_t *bed_filter,
                            const view_output_policy_t *output_policy,
                            const isam_tag_set_t *shared_tags,
                            const isam_tag_set_t *dropped_tags,
                            int count_only, const bam1_t *stream_payload) {
    int selected = passes_record_filter(filter, rec) &&
                   overlaps_bed(bed_filter, rec->tid, rec->pos0, rec->end0);
    if (!selected && !output_policy->unmap_unselected)
        return 0;
    int need_bam = !count_only || view_filter_needs_bam(filter);
    if (!selected && output_policy->unmap_unselected)
        need_bam = 1;
    if (need_bam) {
        const bam1_t *payload = stream_payload;
        if (!stream_payload) {
            store_record_t srec;
            if (!store_lookup_read_id(v->store, rec->store_read_id, &srec))
                die("store_read_id '%lld' missing from store '%s' (run isamtools rehint)",
                    (long long)rec->store_read_id, v->store->store_path);
            int rc;
            if (!v->store->is_bam && v->store->slice_cache)
                rc = store_borrow_payload(v->store, &srec, canon_buf, &payload);
            else {
                rc = payload_cache_fetch(cache, v->store, &srec, canon_buf);
                payload = canon_buf;
            }
            if (rc != 0) die("payload fetch failed for '%s' (rc=%d)", srec.qname, rc);
            store_record_clear(&srec);
        }
        rebuild_from_payload(v, rec, payload,
                             out_buf, shared_tags, dropped_tags);
        if (selected && !passes_bam_filter(filter, hdr, out_buf))
            selected = 0;
        if (!selected && !output_policy->unmap_unselected)
            return 0;
        apply_output_policy(out_buf, output_policy, selected);
    }
    if (count_only) return 1;

    if (sam_write1(out_sam, hdr, out_buf) < 0)
        die("sam_write1 failed");
    return 1;
}

/* Keep the production default while allowing isolated look-ahead benchmarks
 * against the identical frozen decoder/library objects. */
#ifndef ISAM_COORDINATE_WINDOW
#define ISAM_COORDINATE_WINDOW 512
#endif
#if ISAM_COORDINATE_WINDOW < 1 || ISAM_COORDINATE_WINDOW > 65536
#error "ISAM_COORDINATE_WINDOW must be in 1..65536"
#endif
enum { COORDINATE_WINDOW = ISAM_COORDINATE_WINDOW };

static int coordinate_read_record(view_t *v, lmdb_view_iter_t *iterator,
                                  lmdb_view_record_t *rec, store_record_t *store,
                                  view_filter_t *filter,
                                  const view_output_policy_t *policy) {
    if (!lmdb_view_iter_next(iterator, rec)) return 0;
    if (passes_record_filter(filter, rec) || policy->unmap_unselected)
        if (!store_lookup_read_id(v->store, rec->store_read_id, store))
            die("store_read_id '%lld' missing from store '%s'",
                (long long)rec->store_read_id, v->store->store_path);
    return 1;
}

static void coordinate_plan_payloads(view_t *v, lmdb_view_iter_t *iterator,
                                      view_filter_t *filter,
                                      const view_output_policy_t *policy) {
    store_payload_plan_begin(v->store);
    lmdb_view_record_t record;
    store_record_t store = {0};
    while (coordinate_read_record(v, iterator, &record, &store, filter, policy)) {
        if (store.qname) store_payload_plan_add(v->store, &store);
        lmdb_view_record_clear(&record);
        store_record_clear(&store);
    }
    store_payload_plan_finish(v->store);
}

/* Advance the scheduling cursor once per metadata record, not by rescanning
 * the entire resident lookahead on every output record. */
static void coordinate_prefetch_window(store_t *st, store_record_t *stores,
                                       size_t head, size_t count, size_t *cursor,
                                       size_t capacity) {
    while (*cursor < count) {
        size_t start = (head + *cursor) % capacity;
        size_t span = capacity - start;
        if (span > count - *cursor) span = count - *cursor;
        size_t examined = store_prefetch_payloads(st, stores + start, span);
        *cursor += examined;
        if (examined < span) break;
    }
}

/* Rolling coordinate metadata ring. Decode workers stay active while current
 * records are rebuilt/written; only a demanded container can block output. */
static int64_t emit_coordinate_pipeline(samFile *out, view_t *v,
                                       lmdb_view_iter_t *iterator,
                                       bam1_t *scratch, bam1_t *outb,
                                       view_filter_t *filter,
                                       const view_output_policy_t *policy,
                                       const isam_tag_set_t *shared_tags,
                                       const isam_tag_set_t *dropped_tags,
                                       size_t capacity) {
    lmdb_view_record_t *records = xcalloc(capacity, sizeof *records);
    store_record_t *stores = xcalloc(capacity, sizeof *stores);
    int64_t emitted = 0;
    size_t head = 0, count = 0, prefetched = 0;
    int exhausted = 0;
    while (count < capacity) {
        if (!coordinate_read_record(v, iterator, &records[count], &stores[count], filter, policy)) {
            exhausted = 1;
            break;
        }
        ++count;
    }
    while (count) {
        const bam1_t *payload = NULL;
        /* A previously resident payload can have been evicted since the
         * cursor passed it. Give a current miss priority over lookahead. */
        if (stores[head].qname && !store_payload_cached(v->store, &stores[head]))
            store_prefetch_payloads(v->store, &stores[head], 1);
        coordinate_prefetch_window(v->store, stores, head, count, &prefetched, capacity);
        if (stores[head].qname) {
            int rc = store_borrow_payload(v->store, &stores[head], scratch, &payload);
            if (rc) die("payload fetch failed for '%s' (rc=%d)", stores[head].qname, rc);
        }
        /* Consuming a pending result frees its ticket/decoder. Refill before
         * reconstruction and compression, without changing cache residency. */
        coordinate_prefetch_window(v->store, stores, head, count, &prefetched, capacity);
        emitted += emit_lmdb_record(out, v->hdr, v, &records[head], NULL,
                                    scratch, outb, filter, NULL, policy,
                                    shared_tags, dropped_tags, 0, payload);
        if (stores[head].qname) store_payload_release(v->store, &stores[head]);
        lmdb_view_record_clear(&records[head]);
        store_record_clear(&stores[head]);
        head = (head + 1) % capacity;
        --count;
        if (prefetched) --prefetched;
        if (!exhausted) {
            size_t tail = (head + count) % capacity;
            if (coordinate_read_record(v, iterator, &records[tail], &stores[tail], filter, policy))
                ++count;
            else
                exhausted = 1;
        }
    }
    free(stores);
    free(records);
    return emitted;
}

/* Match the exact payload locator, not just a primary flag: the canonical
 * index may use a full-length supplementary-only record or upgrade a preceding
 * supplementary record to a later primary. No retained input is rewritten. */
static int64_t emit_store_order(samFile *out, view_t *v, bam1_t *canon,
                                bam1_t *outb, view_filter_t *filter,
                                const view_output_policy_t *policy,
                                const isam_tag_set_t *shared_tags,
                                const isam_tag_set_t *dropped_tags,
                                int count_only) {
    store_t *st = v->store;
    if (!st->is_bam && (!st->cram_coff || st->cram_nc == 0))
        die("store order requires a canonical CRAM container table; rebuild the store index");
    int64_t raw = 0, ci = 0, base = 0, visited = 0, emitted = 0;
    for (;;) {
        int64_t voffset = st->is_bam ? bgzf_tell(st->sam->fp.bgzf) : -1;
        int rc = sam_read1(st->sam, st->hdr, canon);
        if (rc == -1) break;
        if (rc < -1) die("error reading canonical store '%s' at record %lld",
                         st->store_path, (long long)(raw + 1));
        if (!st->is_bam) {
            while (ci < st->cram_nc && raw >= base + st->cram_nrec[ci])
                base += st->cram_nrec[ci++];
            if (ci >= st->cram_nc)
                die("canonical store has extra records beyond its container table");
        }
        ++raw;
        if (canon->core.flag & BAM_FSECONDARY) continue;
        store_record_t srec;
        if (!store_lookup_qname(st, bam_get_qname(canon), &srec)) {
            if (!(canon->core.flag & BAM_FSUPPLEMENTARY))
                die("canonical primary '%s' is missing from its store index", bam_get_qname(canon));
            continue;
        }
        int matches;
        if (st->is_bam) {
            if (!srec.has_voffset) die("store order requires BAM payload locators; rebuild the store index");
            matches = srec.bam_voffset == voffset;
        } else {
            if (!srec.has_cram_loc) die("store order requires CRAM payload locators; rebuild the store index");
            matches = srec.cram_coffset == st->cram_coff[ci] &&
                      srec.cram_within == raw - 1 - base;
        }
        if (matches) {
            lmdb_view_iter_t *it = lmdb_view_iter_open_read_id(v->lmdb, srec.read_id);
            lmdb_view_record_t rec;
            while (lmdb_view_iter_next(it, &rec)) {
                ++visited;
                emitted += emit_lmdb_record(out, v->hdr, v, &rec, NULL,
                                             NULL, outb, filter, NULL, policy,
                                             shared_tags, dropped_tags,
                                             count_only, canon);
                lmdb_view_record_clear(&rec);
            }
            lmdb_view_iter_close(it);
        }
        store_record_clear(&srec);
    }
    if (!st->is_bam) {
        int64_t expected_raw = 0;
        for (int64_t i = 0; i < st->cram_nc; ++i) expected_raw += st->cram_nrec[i];
        if (raw != expected_raw) die("canonical CRAM record count differs from its container table");
    }
    if (visited != lmdb_view_record_count(v->lmdb))
        die("store order reconstructed %lld view records; expected %lld (stale or incomplete store index)",
            (long long)visited, (long long)lmdb_view_record_count(v->lmdb));
    fprintf(stderr, "isamtools view: store scan %lld canonical records, %lld view records\n",
            (long long)raw, (long long)visited);
    return emitted;
}

/* HTSlib's local REF_PATH route mmaps raw uppercase reference files rather
 * than rereading and normalizing FASTA on every multi-reference cache miss.
 * Require complete coverage of the canonical header, so missing entries cannot
 * silently fall back to the CRAM's original UR path or an online resolver.
 * Files must have been checksum-verified when this directory was populated. */
static void configure_reference_path(store_t *st, const char *directory) {
    struct stat sb;
    if (!*directory || strchr(directory, ':') || strchr(directory, '%') ||
        stat(directory, &sb) != 0 || !S_ISDIR(sb.st_mode))
        die("reference-path '%s' must be an existing directory without ':' or '%%'", directory);
    kstring_t md5 = {0, 0, NULL}, filename = {0, 0, NULL};
    for (int tid = 0; tid < sam_hdr_nref(st->hdr); ++tid) {
        const char *name = sam_hdr_tid2name(st->hdr, tid);
        if (sam_hdr_find_tag_id(st->hdr, "SQ", "SN", name, "M5", &md5) < 0 || md5.l != 32)
            die("canonical reference '%s' lacks a valid M5 tag", name);
        for (size_t i = 0; i < md5.l; ++i) {
            if (!((md5.s[i] >= '0' && md5.s[i] <= '9') ||
                  (md5.s[i] >= 'a' && md5.s[i] <= 'f') ||
                  (md5.s[i] >= 'A' && md5.s[i] <= 'F')))
                die("canonical reference '%s' has a non-hex M5 tag", name);
        }
        filename.l = 0;
        if (ksprintf(&filename, "%s/%s", directory, md5.s) < 0)
            die("cannot allocate reference-path filename");
        if (stat(filename.s, &sb) != 0 || !S_ISREG(sb.st_mode) ||
            sb.st_size != sam_hdr_tid2len(st->hdr, tid))
            die("reference-path file '%s' is missing or has the wrong length for '%s'", filename.s, name);
    }
    filename.l = 0;
    if (ksprintf(&filename, "%s/%%s", directory) < 0 ||
        setenv("REF_PATH", filename.s, 1) != 0 || setenv("REF_CACHE", "", 1) != 0)
        die("cannot configure local CRAM REF_PATH");
    free(md5.s);
    free(filename.s);
}

int cmd_view(int argc, char **argv) {
    const char *store_override = NULL;
    const char *store_idx_override = NULL;
    const char *reference_path = NULL;
    const char *reference_m5_directory = NULL;
    const char *prepare_reference_directory = NULL;
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
    int decode_threads = 0;
    int decode_lookahead = COORDINATE_WINDOW, lookahead_set = 0;
    int cache_selected = 0;
    int write_index = 0;
    const char *order = "coordinate";
    view_filter_t filter = { 0 };
    view_output_policy_t output_policy = { 0 };
    isam_tag_set_init(&output_policy.remove_tags);
    isam_tag_set_init(&output_policy.keep_tags);

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
        {"reference-path", required_argument, 0, 1017},
        {"prepare-reference", required_argument, 0, 1019},
        {"threads",     required_argument, 0, '@'},
        {"decode-threads", required_argument, 0, 1018},
        {"decode-lookahead", required_argument, 0, 1020},
        {"cache-selected", no_argument, 0, 1021},
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
                            "1bBCcd:D:e:f:F:G:hHL:l:Mm:N:O:o:Ppq:R:r:Ss:t:T:uU:Xx:z:@:?",
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
        case 'p': output_policy.unmap_unselected = 1; break;
        case 'q': filter.min_mapq = atoi(optarg); break;
        case 'R': str_set_load_file(&filter.read_groups, optarg, "read-group"); break;
        case 'r': str_set_add(&filter.read_groups, optarg); break;
        case 'S': break;
        case 's': parse_subsample_s(&filter, optarg); break;
        case 't': break;
        case 'T': reference_path = optarg; break;
        case 'u': out_fmt = VIEW_OUT_BAM; compression_level = 0; break;
        case 'U': unsupported_option("-U/--output-unselected"); break;
        case 'X': unsupported_option("-X/--customized-index"); break;
        case 'x':
            if (optarg[0] == '^') {
                parse_tag_list_into(&output_policy.keep_tags, optarg + 1);
                output_policy.keep_tags_active = 1;
            } else {
                parse_tag_list_into(&output_policy.remove_tags, optarg);
            }
            break;
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
        case 1009:
            filter.subsample_fraction = parse_fraction(optarg, "subsample fraction");
            filter.subsample_active = 1;
            break;
        case 1010:
            filter.subsample_seed = (uint64_t)parse_nonnegative_int(optarg, "subsample seed");
            break;
        case 1011: unsupported_option("--add-flags"); break;
        case 1012: unsupported_option("--remove-flags"); break;
        case 1013:
            parse_tag_list_into(&output_policy.keep_tags, optarg);
            output_policy.keep_tags_active = 1;
            break;
        case 1014: write_index = 1; break;
        case 1015: break;
        case 1016: unsupported_option("--input-fmt-option"); break;
        case 1017: reference_m5_directory = optarg; break;
        case 1019: prepare_reference_directory = optarg; break;
        case 1020:
            decode_lookahead = parse_nonnegative_int(optarg, "decode look-ahead");
            if (decode_lookahead < 1 || decode_lookahead > 65536)
                die("decode look-ahead must be in 1..65536");
            lookahead_set = 1;
            break;
        case 1021: cache_selected = 1; break;
        case 1018: {
            int64_t workers = parse_nonnegative_int(optarg, "decode threads");
            if (workers < 1 || workers > 64)
                die("decode threads must be in 1..64");
            decode_threads = (int)workers;
            break;
        }
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *isam_path = argv[optind++];
    int n_regions = argc - optind;
    char **regions = argv + optind;
    const char *region = n_regions == 1 ? regions[0] : NULL;
    if (strcasecmp(order, "coordinate") && strcasecmp(order, "name") && strcasecmp(order, "store"))
        die("unknown view order '%s' (expected coordinate, name, or store)", order);
    int store_order = strcasecmp(order, "store") == 0;
    if (store_order && (n_regions || bed_path || subset || fetch_pairs))
        die("store order is a full scan and cannot be combined with regions, BED, subsets or fetch-pairs");
    if (decode_threads && (strcasecmp(order, "coordinate") || n_regions > 1 ||
                           bed_path || subset || fetch_pairs || count_only ||
                           header_only || output_policy.unmap_unselected))
        die("decode threads require coordinate output, whole view or one region, without BED/subset/fetch-pairs/count/unmap");
    if ((lookahead_set || cache_selected) && !decode_threads)
        die("--decode-lookahead and --cache-selected require --decode-threads");
    if (write_index && (out_fmt != VIEW_OUT_BAM || !strcmp(output_path, "-") ||
                        strcasecmp(order, "coordinate") || n_regions > 1 ||
                        bed_path || subset || fetch_pairs || count_only ||
                        header_only || output_policy.unmap_unselected))
        die("--write-index requires a named coordinate BAM, whole view or one region, without BED/subset/fetch-pairs/count/header/unmap");
    if (write_index) {
        struct stat output_stat;
        if (stat(output_path, &output_stat) == 0 && !S_ISREG(output_stat.st_mode))
            die("--write-index requires a regular BAM output file");
    }
    if (out_fmt == VIEW_OUT_CRAM && !reference_path)
        die("CRAM output requires -T/--reference");
    if (reference_path && reference_m5_directory)
        die("choose either -T/--reference or --reference-path, not both");
    if (prepare_reference_directory && (!reference_path || reference_m5_directory ||
                                      count_only || header_only || out_fmt == VIEW_OUT_CRAM))
        die("--prepare-reference requires -T and SAM/BAM record output, without --reference-path");
    if (out_fmt == VIEW_OUT_SAM && compression_level >= 0)
        die("compression options require BAM or CRAM output");
    if (fetch_pairs && n_regions == 0 && !bed_path)
        die("-P/--fetch-pairs requires a region or -L/--targets-file");
    if (fetch_pairs && output_policy.unmap_unselected)
        die("-P/--fetch-pairs cannot be combined with -p/--unmap");

    view_t *v = view_open_read(isam_path);
    if (count_only && !header_only && n_regions == 0 && !bed_path &&
        !subset && !qname_path && !fetch_pairs && view_filter_is_noop(&filter)) {
        int64_t count = lmdb_view_record_count(v->lmdb);
        write_count(output_path, count);
        view_filter_destroy(&filter);
        view_output_policy_destroy(&output_policy);
        view_close(v);
        return 0;
    }
    view_link_store(v, store_override, store_idx_override);
    if (prepare_reference_directory) {
        if (v->store->is_bam) die("--prepare-reference requires a canonical CRAM store");
        if (isam_reference_prepare(v->store->hdr, reference_path, prepare_reference_directory) < 0)
            die("raw-reference preparation failed");
        reference_m5_directory = prepare_reference_directory;
        reference_path = NULL;
    }
    if (!v->store->is_bam && !reference_path && !reference_m5_directory)
        die("canonical store '%s' is CRAM; -T/--reference or --reference-path is required to decode read payload",
            v->store->store_path);
    if (reference_m5_directory) {
        if (v->store->is_bam) die("--reference-path requires a canonical CRAM store");
        configure_reference_path(v->store, reference_m5_directory);
    }
    if (!v->store->is_bam && reference_path && !store_order && !count_only && !header_only)
        warnx_("FASTA coordinate decoding can be slower. Prepare raw references with "
               "'isamtools ref-index --cram STORE -T FASTA -o DIR', then use "
               "'--reference-path DIR'; or add '--prepare-reference DIR' for one-step preparation.");
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
    /* Container-locator reads seek after each cache miss. HTSlib's threaded
     * CRAM reader prefetches additional slices, then cram_seek drains and
     * discards that work. Decode only the requested container synchronously;
     * retain parallel I/O for BAM/legacy readers and output compression. */
    if (n_threads > 1 && (store_order || v->store->is_bam || !v->store->slice_cache) &&
        hts_set_threads(v->store->sam, n_threads) != 0)
        die("cannot set store threads to %d", n_threads);
    if (decode_threads)
        store_set_decode_threads(v->store, decode_threads, reference_path);
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

    if (store_order) {
        if (sam_hdr_update_hd(v->hdr, "SO", "unsorted") != 0 ||
            sam_hdr_remove_tag_hd(v->hdr, "SS") < 0)
            die("cannot mark store-order output header as unsorted");
    }

    if (out && (print_header || out_fmt != VIEW_OUT_SAM)) {
        if (sam_hdr_write(out, v->hdr) < 0) die("sam_hdr_write failed");
    }
    kstring_t index_path = {0, 0, NULL};
    if (write_index) {
        if (ksprintf(&index_path, "%s.bai", output_path) < 0 ||
            sam_idx_init(out, v->hdr, 0, index_path.s) < 0)
            die("cannot initialize BAM output index '%s'", output_path);
    }
    if (header_only) {
        bed_filter_destroy(&bed_filter);
        view_filter_destroy(&filter);
        view_output_policy_destroy(&output_policy);
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
            if (output_policy.unmap_unselected)
                die("-p/--unmap cannot be combined with '*' unmapped region selection");
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
            if (output_policy.unmap_unselected) {
                if (use_bed_filter) {
                    bed_filter_t intersected = { 0 };
                    for (size_t i = 0; i < bed_filter.n; ++i) {
                        if (bed_filter.items[i].tid != tid)
                            continue;
                        int64_t ibeg = bed_filter.items[i].beg > beg
                            ? bed_filter.items[i].beg : beg;
                        int64_t iend = bed_filter.items[i].end < end
                            ? bed_filter.items[i].end : end;
                        bed_filter_add(&intersected, tid, ibeg, iend);
                    }
                    bed_filter_destroy(&bed_filter);
                    bed_filter = intersected;
                } else {
                    bed_filter_add(&bed_filter, tid, beg, end);
                    use_bed_filter = 1;
                }
                region = NULL;
            } else {
                iter_kind = LMDB_VIEW_ITER_REGION;
            }
        }
    }

    if (store_order) {
        n_emit = emit_store_order(out, v, canon, outb, &filter, &output_policy,
                                  &shared_tags, &dropped_tags, count_only);
    } else if (!output_policy.unmap_unselected &&
        have_subset && !subset_from_store && (use_bed_filter || region)) {
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
                                            &output_policy,
                                            &shared_tags, &dropped_tags,
                                            fetch_pairs ? 1 : count_only, NULL);
            if (fetch_pairs) {
                if (selected) i64_set_add(&pair_reads, rec.store_read_id);
            } else {
                n_emit += selected;
            }
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
    } else if (use_bed_filter && !output_policy.unmap_unselected) {
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
                                                    &output_policy,
                                                    &shared_tags, &dropped_tags,
                                                    fetch_pairs ? 1 : count_only, NULL);
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
        if (count_only && !region && !use_bed_filter &&
            (!have_subset || subset_from_store)) {
            iter_kind = LMDB_VIEW_ITER_RECORD;
        } else if (!region && have_subset && !subset_from_store) {
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
        if (decode_threads) {
            if (iter_kind != LMDB_VIEW_ITER_COORD && iter_kind != LMDB_VIEW_ITER_REGION &&
                iter_kind != LMDB_VIEW_ITER_UNMAPPED)
                die("decode threads require a coordinate-indexed view");
            if (cache_selected) {
                coordinate_plan_payloads(v, it, &filter, &output_policy);
                lmdb_view_iter_close(it);
                it = lmdb_view_iter_open(v->lmdb, iter_kind, tid, beg, end);
            }
            n_emit += emit_coordinate_pipeline(out, v, it, canon, outb, &filter,
                                              &output_policy, &shared_tags, &dropped_tags,
                                              (size_t)decode_lookahead);
        } else while (lmdb_view_iter_next(it, &rec)) {
            if (have_subset && subset_from_store &&
                !subset_contains(v, subset_from_store, subset_id,
                                 rec.store_read_id)) {
                lmdb_view_record_clear(&rec);
                continue;
            }
            int selected = emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                            canon, outb, &filter,
                                            use_bed_filter ? &bed_filter : NULL,
                                            &output_policy,
                                            &shared_tags, &dropped_tags,
                                            fetch_pairs ? 1 : count_only, NULL);
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
                                               &output_policy,
                                               &shared_tags, &dropped_tags,
                                               count_only, NULL);
                lmdb_view_record_clear(&rec);
            }
            lmdb_view_iter_close(it);
        }
        seen_aln_ids_destroy(&emitted);
    }

    if (count_only) {
        write_count(output_path, n_emit);
    }
    i64_set_destroy(&pair_reads);
    bed_filter_destroy(&bed_filter);
    view_filter_destroy(&filter);
    view_output_policy_destroy(&output_policy);
    payload_cache_destroy(cache);
    isam_tag_set_destroy(&shared_tags);
    isam_tag_set_destroy(&dropped_tags);
    bam_destroy1(canon);
    bam_destroy1(outb);
    if (write_index && sam_idx_save(out) < 0)
        die("cannot save BAM output index '%s'", index_path.s);
    if (out && sam_close(out) < 0) die("cannot close output '%s'", output_path);
    free(index_path.s);
    view_close(v);
    return 0;
}
