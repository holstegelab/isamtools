#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <errno.h>
#include <strings.h>
#include <sys/types.h>

static void usage(void) {
    fputs(
        "Usage: isamtools view [options] <view.isam> [region]\n"
        "\n"
        "Emit reconstructed SAM/BAM by combining alignment-view records with\n"
        "the canonical store payload.\n"
        "\n"
        "Options:\n"
        "  -b               Output BAM (default: SAM).\n"
        "  -O, --output-fmt FMT  Output format: sam, bam, or cram.\n"
        "  -h, --header     Include @ header lines (default).\n"
        "  -H, --header-only  Print only the header and exit.\n"
        "      --no-header  Suppress the header.\n"
        "      --store P    Override store path.\n"
        "      --store-index P  Override store index path.\n"
        "  -T, --reference P Reference FASTA for CRAM payload decoding.\n"
        "  -@, --threads N  Worker threads for htslib I/O.\n"
        "  -L, --targets-file P  BED file of regions to emit.\n"
        "      --order ORD  full-scan order: coordinate (default) or name.\n"
        "      --subset S   Only emit reads in named view subset.\n"
        "  -f INT           Require all FLAG bits in INT.\n"
        "  -F INT           Exclude records with any FLAG bits in INT.\n"
        "      --rf INT     Require at least one FLAG bit in INT.\n"
        "  -G INT           Exclude records with all FLAG bits in INT.\n"
        "  -q INT           Minimum MAPQ.\n"
        "  -c, --count      Count records instead of emitting them.\n"
        "  -o, --output P   Output file (default stdout).\n"
        "  --help           Show this help.\n",
        stderr);
}

typedef struct view_filter {
    int require_all;
    int exclude_any;
    int require_any;
    int exclude_all;
    int min_mapq;
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

static void bed_filter_add(bed_filter_t *bed, int tid,
                           int64_t beg, int64_t end) {
    if (end <= beg) return;
    if (bed->n == bed->cap) {
        bed->cap = bed->cap ? bed->cap * 2 : 128;
        bed->items = xrealloc(bed->items, bed->cap * sizeof *bed->items);
    }
    bed->items[bed->n++] = (bed_interval_t){ tid, beg, end };
}

static bed_filter_t bed_filter_load(const char *path, sam_hdr_t *hdr) {
    bed_filter_t bed = { 0 };
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
        bed_filter_add(&bed, tid, beg, end);
    }
    free(line);
    fclose(fp);
    qsort(bed.items, bed.n, sizeof *bed.items, bed_interval_cmp);
    return bed;
}

static void bed_filter_destroy(bed_filter_t *bed) {
    free(bed->items);
    memset(bed, 0, sizeof *bed);
}

static int overlaps_bed(const bed_filter_t *bed,
                        int tid, int64_t beg, int64_t end) {
    if (!bed) return 1;
    if (bed->n == 0) return 0;
    if (tid < 0 || beg < 0 || end <= beg) return 0;
    for (size_t i = 0; i < bed->n; ++i) {
        const bed_interval_t *iv = &bed->items[i];
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

static int parse_flag_arg(const char *s) {
    char *end = NULL;
    long v = strtol(s, &end, 0);
    if (!s[0] || (end && *end) || v < 0 || v > 0xffff)
        die("invalid flag value '%s'", s);
    return (int)v;
}

static int passes_filter(const view_filter_t *f, int flag, int mapq) {
    if (f->require_all && (flag & f->require_all) != f->require_all) return 0;
    if (f->exclude_any && (flag & f->exclude_any) != 0) return 0;
    if (f->require_any && (flag & f->require_any) == 0) return 0;
    if (f->exclude_all && (flag & f->exclude_all) == f->exclude_all) return 0;
    if (mapq < f->min_mapq) return 0;
    return 1;
}

static int emit_lmdb_record(samFile *out_sam, sam_hdr_t *hdr,
                            view_t *v, const lmdb_view_record_t *rec,
                            payload_cache_t *cache,
                            bam1_t *canon_buf, bam1_t *out_buf,
                            const view_filter_t *filter,
                            const bed_filter_t *bed,
                            int count_only) {
    if (!passes_filter(filter, rec->flag, rec->mapq)) return 0;
    if (!overlaps_bed(bed, rec->tid, rec->pos0, rec->end0)) return 0;
    if (count_only) return 1;

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
    store_record_clear(&srec);

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
    int mode_bam = 0;
    int mode_cram = 0;
    int print_header = 1;
    int header_only  = 0;
    int count_only = 0;
    int n_threads = 0;
    const char *order = "coordinate";
    view_filter_t filter = { 0, 0, 0, 0, 0 };

    static struct option longopts[] = {
        {"header",      no_argument, 0, 'h'},
        {"header-only", no_argument, 0, 'H'},
        {"no-header",   no_argument, 0, 1000},
        {"store",       required_argument, 0, 1001},
        {"store-index", required_argument, 0, 1002},
        {"order",       required_argument, 0, 1003},
        {"output",      required_argument, 0, 'o'},
        {"output-fmt",  required_argument, 0, 'O'},
        {"help",        no_argument, 0, 1004},
        {"rf",          required_argument, 0, 1005},
        {"subset",      required_argument, 0, 1006},
        {"targets-file", required_argument, 0, 'L'},
        {"count",       no_argument,       0, 'c'},
        {"reference",   required_argument, 0, 'T'},
        {"threads",     required_argument, 0, '@'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "bcf:F:G:hHL:O:o:q:T:@:", longopts, NULL)) != -1) {
        switch (c) {
        case 'b': mode_bam = 1; mode_cram = 0; break;
        case 'c': count_only = 1; print_header = 0; break;
        case 'f': filter.require_all = parse_flag_arg(optarg); break;
        case 'F': filter.exclude_any = parse_flag_arg(optarg); break;
        case 'G': filter.exclude_all = parse_flag_arg(optarg); break;
        case 'h': print_header = 1; break;
        case 'H': header_only = 1; print_header = 1; break;
        case 'L': bed_path = optarg; break;
        case 'O':
            if (strcasecmp(optarg, "sam") == 0) {
                mode_bam = 0;
                mode_cram = 0;
            } else if (strcasecmp(optarg, "bam") == 0) {
                mode_bam = 1;
                mode_cram = 0;
            } else if (strcasecmp(optarg, "cram") == 0) {
                mode_bam = 0;
                mode_cram = 1;
            } else {
                die("unknown output format '%s' (expected sam, bam, or cram)", optarg);
            }
            break;
        case 'q': filter.min_mapq = atoi(optarg); break;
        case 'T': reference_path = optarg; break;
        case '@': n_threads = atoi(optarg); break;
        case 1000: print_header = 0; break;
        case 1001: store_override = optarg; break;
        case 1002: store_idx_override = optarg; break;
        case 1003: order = optarg; break;
        case 'o': output_path = optarg; break;
        case 1004: usage(); return 0;
        case 1005: filter.require_any = parse_flag_arg(optarg); break;
        case 1006: subset = optarg; break;
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *isam_path = argv[optind++];
    const char *region    = (optind < argc) ? argv[optind] : NULL;

    view_t *v = view_open_read(isam_path);
    view_link_store(v, store_override, store_idx_override);
    bed_filter_t bed = { 0 };
    if (bed_path) bed = bed_filter_load(bed_path, v->hdr);
    if (reference_path &&
        hts_set_fai_filename(v->store->sam, reference_path) != 0)
        die("cannot set reference FASTA '%s'", reference_path);
    if (n_threads > 1 &&
        hts_set_threads(v->store->sam, n_threads) != 0)
        die("cannot set store threads to %d", n_threads);
    int have_subset = 0;
    int64_t subset_id = 0;
    if (subset) {
        if (!lmdb_subset_lookup(v->lmdb, subset, &subset_id))
            die("subset '%s' is empty or absent in view", subset);
        have_subset = 1;
    }

    /* Open output. */
    char mode_str[8];
    snprintf(mode_str, sizeof mode_str, "w%s", mode_cram ? "c" : mode_bam ? "b" : "");
    samFile *out = sam_open(output_path, mode_str);
    if (!out) die("cannot open output '%s'", output_path);
    if (n_threads > 1 && hts_set_threads(out, n_threads) != 0)
        die("cannot set output threads to %d", n_threads);

    if (print_header) {
        if (sam_hdr_write(out, v->hdr) < 0) die("sam_hdr_write failed");
    }
    if (header_only) {
        bed_filter_destroy(&bed);
        sam_close(out);
        view_close(v);
        return 0;
    }

    bam1_t *canon = bam_init1();
    bam1_t *outb  = bam_init1();
    payload_cache_t *cache = payload_cache_create(64 * 1024 * 1024);
    int64_t n_emit = 0;

    lmdb_view_iter_kind_t iter_kind = LMDB_VIEW_ITER_IMPORT;
    int tid = -1;
    int64_t beg = 0, end = 0;
    if (region) {
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

    if (have_subset && (bed_path || region)) {
        lmdb_view_iter_t *it = lmdb_view_iter_open_subset(v->lmdb, subset_id);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            if (region && !(rec.tid == tid && rec.pos0 < end && rec.end0 > beg)) {
                lmdb_view_record_clear(&rec);
                continue;
            }
            n_emit += emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                       canon, outb, &filter,
                                       bed_path ? &bed : NULL, count_only);
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
    } else if (bed_path) {
        int64_t has_iv = 0;
        view_meta_get_i64(v, "indexed_interval", &has_iv);
        if (!has_iv)
            die("alignment view has no LMDB interval index; "
                "run `isamtools sort --interval %s`", isam_path);
        seen_aln_ids_t seen = { 0 };
        for (size_t i = 0; i < bed.n; ++i) {
            int iv_tid = bed.items[i].tid;
            int64_t iv_beg = bed.items[i].beg;
            int64_t iv_end = bed.items[i].end;
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
                if (!seen_aln_id(&seen, rec.aln_id))
                    n_emit += emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                               canon, outb, &filter, &bed,
                                               count_only);
                lmdb_view_record_clear(&rec);
            }
            lmdb_view_iter_close(it);
        }
        seen_aln_ids_destroy(&seen);
    } else {
        if (!region && have_subset) {
            iter_kind = LMDB_VIEW_ITER_SUBSET;
        } else if (!region && strcasecmp(order, "name") == 0) {
            iter_kind = LMDB_VIEW_ITER_NAME;
        } else if (!region) {
            int64_t has_coord = 0;
            view_meta_get_i64(v, "indexed_coord", &has_coord);
            iter_kind = has_coord ? LMDB_VIEW_ITER_COORD : LMDB_VIEW_ITER_IMPORT;
        }

        lmdb_view_iter_t *it = have_subset && !region
            ? lmdb_view_iter_open_subset(v->lmdb, subset_id)
            : lmdb_view_iter_open(v->lmdb, iter_kind, tid, beg, end);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            n_emit += emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                       canon, outb, &filter,
                                       bed_path ? &bed : NULL, count_only);
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);
    }

    if (count_only) fprintf(stdout, "%lld\n", (long long)n_emit);
    bed_filter_destroy(&bed);
    payload_cache_destroy(cache);
    bam_destroy1(canon);
    bam_destroy1(outb);
    sam_close(out);
    view_close(v);
    return 0;
}
