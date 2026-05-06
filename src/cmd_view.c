#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <strings.h>

static void usage(void) {
    fputs(
        "Usage: isamtools view [options] <view.isam> [region]\n"
        "\n"
        "Emit reconstructed SAM/BAM by combining alignment-view records with\n"
        "the canonical store payload.\n"
        "\n"
        "Options:\n"
        "  -b               Output BAM (default: SAM).\n"
        "  -h, --header     Include @ header lines (default).\n"
        "  -H, --header-only  Print only the header and exit.\n"
        "      --no-header  Suppress the header.\n"
        "      --store P    Override store path.\n"
        "      --store-index P  Override store index path.\n"
        "  -T, --reference P Reference FASTA for CRAM payload decoding.\n"
        "  -@, --threads N  Worker threads for htslib I/O.\n"
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

static void prepare_subset(sqlite3 *db, const char *subset) {
    if (!subset) return;
    db_exec(db, "DROP TABLE IF EXISTS temp.selected_subset;");
    db_exec(db,
        "CREATE TEMP TABLE selected_subset("
        "  store_read_id INTEGER PRIMARY KEY"
        ");");
    sqlite3_stmt *s = db_prep(db,
        "INSERT INTO selected_subset(store_read_id) "
        "SELECT sr.store_read_id "
        "FROM subset s JOIN subset_read sr ON sr.subset_id=s.subset_id "
        "WHERE s.name=?;");
    sqlite3_bind_text(s, 1, subset, -1, SQLITE_STATIC);
    db_step_done(s);
    sqlite3_finalize(s);
    int64_t n;
    sqlite3_stmt *c = db_prep(db, "SELECT COUNT(*) FROM selected_subset;");
    if (sqlite3_step(c) != SQLITE_ROW)
        die("subset count failed: %s", sqlite3_errmsg(db));
    n = sqlite3_column_int64(c, 0);
    sqlite3_finalize(c);
    if (n == 0) die("subset '%s' is empty or absent in view", subset);
}

/* Emit one row's reconstructed bam record. cigar/aux blobs are columns. */
static int emit_row(samFile *out_sam, sam_hdr_t *hdr,
                    view_t *v, sqlite3_stmt *q,
                    payload_cache_t *cache,
                    bam1_t *canon_buf, bam1_t *out_buf,
                    const view_filter_t *filter,
                    int count_only) {
    /* Column order: store_read_id, qname, flag, tid, pos0, mapq, n_cigar, cigar,
     *               mate_tid, mate_pos0, tlen, aux */
    int64_t store_read_id = sqlite3_column_int64(q, 0);
    int flag    = sqlite3_column_int(q, 2);
    int tid     = sqlite3_column_type(q,3)==SQLITE_NULL ? -1 : sqlite3_column_int(q,3);
    int64_t pos = sqlite3_column_type(q,4)==SQLITE_NULL ? -1 : sqlite3_column_int64(q,4);
    int mapq    = sqlite3_column_type(q,5)==SQLITE_NULL ?  0 : sqlite3_column_int(q,5);
    int n_cigar = sqlite3_column_type(q,6)==SQLITE_NULL ?  0 : sqlite3_column_int(q,6);
    const uint32_t *cigar = (const uint32_t *)sqlite3_column_blob(q, 7);
    int mtid    = sqlite3_column_type(q,8)==SQLITE_NULL ? -1 : sqlite3_column_int(q,8);
    int64_t mpos= sqlite3_column_type(q,9)==SQLITE_NULL ? -1 : sqlite3_column_int64(q,9);
    int64_t tlen= sqlite3_column_int64(q, 10);
    const uint8_t *aux = (const uint8_t *)sqlite3_column_blob(q, 11);
    int aux_len = sqlite3_column_bytes(q, 11);

    if (!passes_filter(filter, flag, mapq)) return 0;
    if (count_only) return 1;

    store_record_t srec;
    if (!store_lookup_read_id(v->store, store_read_id, &srec))
        die("store_read_id '%lld' missing from store '%s' (run isamtools rehint)",
            (long long)store_read_id, v->store->store_path);
    int rc = payload_cache_fetch(cache, v->store, &srec, canon_buf);
    if (rc != 0) die("payload fetch failed for '%s' (rc=%d)", srec.qname, rc);

    isam_rebuild_bam(out_buf, canon_buf,
                     flag, tid, pos, mapq,
                     n_cigar, cigar, mtid, mpos, tlen, aux, aux_len);
    store_record_clear(&srec);

    if (sam_write1(out_sam, hdr, out_buf) < 0)
        die("sam_write1 failed");
    return 1;
}

#ifdef ISAM_ENABLE_LMDB
static int emit_lmdb_record(samFile *out_sam, sam_hdr_t *hdr,
                            view_t *v, const lmdb_view_record_t *rec,
                            payload_cache_t *cache,
                            bam1_t *canon_buf, bam1_t *out_buf,
                            const view_filter_t *filter,
                            int count_only) {
    if (!passes_filter(filter, rec->flag, rec->mapq)) return 0;
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
#endif

int cmd_view(int argc, char **argv) {
    const char *store_override = NULL;
    const char *store_idx_override = NULL;
    const char *reference_path = NULL;
    const char *output_path = "-";
    const char *subset = NULL;
    int mode_bam = 0;
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
        {"help",        no_argument, 0, 1004},
        {"rf",          required_argument, 0, 1005},
        {"subset",      required_argument, 0, 1006},
        {"count",       no_argument,       0, 'c'},
        {"reference",   required_argument, 0, 'T'},
        {"threads",     required_argument, 0, '@'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "bcf:F:G:hHo:q:T:@:", longopts, NULL)) != -1) {
        switch (c) {
        case 'b': mode_bam = 1; break;
        case 'c': count_only = 1; print_header = 0; break;
        case 'f': filter.require_all = parse_flag_arg(optarg); break;
        case 'F': filter.exclude_any = parse_flag_arg(optarg); break;
        case 'G': filter.exclude_all = parse_flag_arg(optarg); break;
        case 'h': print_header = 1; break;
        case 'H': header_only = 1; print_header = 1; break;
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
    if (reference_path &&
        hts_set_fai_filename(v->store->sam, reference_path) != 0)
        die("cannot set reference FASTA '%s'", reference_path);
    if (n_threads > 1 &&
        hts_set_threads(v->store->sam, n_threads) != 0)
        die("cannot set store threads to %d", n_threads);
    if (v->backend_kind == VIEW_BACKEND_SQLITE)
        prepare_subset(v->db, subset);
    else if (subset)
        die("LMDB alignment views do not yet support named subsets");

    /* Open output. */
    char mode_str[8];
    snprintf(mode_str, sizeof mode_str, "w%s", mode_bam ? "b" : "");
    samFile *out = sam_open(output_path, mode_str);
    if (!out) die("cannot open output '%s'", output_path);
    if (n_threads > 1 && hts_set_threads(out, n_threads) != 0)
        die("cannot set output threads to %d", n_threads);

    if (print_header) {
        if (sam_hdr_write(out, v->hdr) < 0) die("sam_hdr_write failed");
    }
    if (header_only) { sam_close(out); view_close(v); return 0; }

    bam1_t *canon = bam_init1();
    bam1_t *outb  = bam_init1();
    payload_cache_t *cache = payload_cache_create(64 * 1024 * 1024);
    int64_t n_emit = 0;

#ifdef ISAM_ENABLE_LMDB
    if (v->backend_kind == VIEW_BACKEND_LMDB) {
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
        } else if (strcasecmp(order, "name") == 0) {
            iter_kind = LMDB_VIEW_ITER_NAME;
        } else {
            int64_t has_coord = 0;
            view_meta_get_i64(v, "indexed_coord", &has_coord);
            iter_kind = has_coord ? LMDB_VIEW_ITER_COORD : LMDB_VIEW_ITER_IMPORT;
        }

        lmdb_view_iter_t *it = lmdb_view_iter_open(v->lmdb, iter_kind,
                                                   tid, beg, end);
        lmdb_view_record_t rec;
        while (lmdb_view_iter_next(it, &rec)) {
            n_emit += emit_lmdb_record(out, v->hdr, v, &rec, cache,
                                       canon, outb, &filter, count_only);
            lmdb_view_record_clear(&rec);
        }
        lmdb_view_iter_close(it);

        if (count_only) fprintf(stdout, "%lld\n", (long long)n_emit);
        payload_cache_destroy(cache);
        bam_destroy1(canon);
        bam_destroy1(outb);
        sam_close(out);
        view_close(v);
        return 0;
    }
#endif

    sqlite3_stmt *q = NULL;
    const char *subset_join = subset
        ? "JOIN selected_subset ss ON ss.store_read_id=r.store_read_id "
        : "";

    if (region) {
        /* Parse region against view header. */
        int tid = -1;
        hts_pos_t beg = 0, end = HTS_POS_MAX;
        const char *tail = sam_parse_region(v->hdr, region, &tid,
                                            &beg, &end, 0);
        if (!tail || tid < 0) die("invalid region: '%s'", region);

        /* Require interval index. */
        int64_t has_iv = 0;
        meta_get_i64(v->db, "indexed_interval", &has_iv);
        if (!has_iv)
            die("alignment view has no R*Tree interval index; "
                "run `isamtools sort --interval %s`", isam_path);

        char sql[1400];
        snprintf(sql, sizeof sql,
            "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
            "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
            "FROM aln_interval i "
            "JOIN aln_id_lookup k ON k.aln_id=i.aln_id "
            "JOIN aln_record    r ON r.tid=k.tid AND r.pos0=k.pos0 "
            "                    AND r.record_order=k.record_order "
            "                    AND r.aln_id=k.aln_id "
            "JOIN read_dict     d ON d.store_read_id=r.store_read_id "
            "%s"
            "WHERE i.tid_min<=? AND i.tid_max>=? "
            "  AND i.pos_min<? AND i.pos_max>? "
            "ORDER BY r.tid, r.pos0, r.record_order, r.aln_id;",
            subset_join);
        q = db_prep(v->db, sql);
        sqlite3_bind_int  (q, 1, tid);
        sqlite3_bind_int  (q, 2, tid);
        sqlite3_bind_int64(q, 3, (int64_t)end);
        sqlite3_bind_int64(q, 4, (int64_t)beg);
    } else {
        int64_t has_coord = 0;
        meta_get_i64(v->db, "indexed_coord", &has_coord);
        if (!has_coord) {
            const char *order_clause =
                strcasecmp(order, "name") == 0
                ? "ORDER BY d.qname_hash, d.qname, r.record_order, r.aln_id"
                : "ORDER BY r.is_unmapped, r.tid, r.pos0, r.record_order, r.aln_id";
            char sql[768];
            snprintf(sql, sizeof sql,
                "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
                "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
                "FROM aln_record_import r "
                "JOIN read_dict d ON d.store_read_id=r.store_read_id "
                "%s%s;",
                subset_join,
                order_clause);
            q = db_prep(v->db, sql);
        } else if (strcasecmp(order, "name") == 0) {
            char sql[1400];
            snprintf(sql, sizeof sql,
                "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
                "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
                "FROM ("
                "  SELECT store_read_id, flag, tid, pos0, mapq, n_cigar, cigar,"
                "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id "
                "  FROM aln_record "
                "  UNION ALL "
                "  SELECT store_read_id, flag, NULL, NULL, mapq, n_cigar, cigar,"
                "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id "
                "  FROM aln_unmapped"
                ") r "
                "JOIN read_dict d ON d.store_read_id=r.store_read_id "
                "%s"
                "ORDER BY d.qname_hash, d.qname, r.record_order, r.aln_id;",
                subset_join);
            q = db_prep(v->db, sql);
        } else {
            char sql[768];
            snprintf(sql, sizeof sql,
                "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
                "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
                "FROM aln_record r "
                "JOIN read_dict d ON d.store_read_id=r.store_read_id "
                "%s"
                "ORDER BY r.tid, r.pos0, r.record_order, r.aln_id;",
                subset_join);
            q = db_prep(v->db, sql);
        }
    }

    while (db_step(q) == SQLITE_ROW)
        n_emit += emit_row(out, v->hdr, v, q, cache, canon, outb,
                           &filter, count_only);
    sqlite3_finalize(q);

    /* No region + no order=name: also emit unmapped records appended,
     * matching samtools-style behavior. */
    int64_t has_coord = 0;
    meta_get_i64(v->db, "indexed_coord", &has_coord);
    if (!region && has_coord && strcasecmp(order, "name") != 0) {
        char sql[900];
        snprintf(sql, sizeof sql,
            "SELECT u.store_read_id, d.qname, u.flag, NULL, NULL, u.mapq, u.n_cigar, u.cigar,"
            "       u.mate_tid, u.mate_pos0, u.tlen, u.aux "
            "FROM aln_unmapped u JOIN read_dict d ON d.store_read_id=u.store_read_id "
            "%s"
            "ORDER BY u.record_order, u.aln_id;",
            subset ? "JOIN selected_subset ss ON ss.store_read_id=u.store_read_id " : "");
        sqlite3_stmt *qu = db_prep(v->db, sql);
        while (db_step(qu) == SQLITE_ROW)
            n_emit += emit_row(out, v->hdr, v, qu, cache, canon, outb,
                               &filter, count_only);
        sqlite3_finalize(qu);
    }

    if (count_only) fprintf(stdout, "%lld\n", (long long)n_emit);
    payload_cache_destroy(cache);
    bam_destroy1(canon);
    bam_destroy1(outb);
    sam_close(out);
    view_close(v);
    return 0;
}
