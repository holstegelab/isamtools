#include "common.h"
#include "schema.h"
#include "backend_lmdb.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>
#include <inttypes.h>

static void usage(void) {
    fputs(
        "Usage: isamtools index [options] <store.bam|store.cram>\n"
        "\n"
        "Build the canonical store index <store>.isam.\n"
        "By default builds the QNAME index (required for `import`).\n"
        "\n"
        "Options:\n"
        "  --name NAME      Logical alias for this store (e.g. chm13).\n"
        "  --backend NAME   Backend: sqlite (default) or lmdb.\n"
        "  -o, --output P   Override default output path (<store>.isam).\n"
        "  -h, --help       Show this help.\n",
        stderr);
}

static int64_t scalar_i64(sqlite3 *db, const char *sql) {
    sqlite3_stmt *s = db_prep(db, sql);
    int rc = sqlite3_step(s);
    if (rc != SQLITE_ROW)
        die("sqlite scalar query returned no row: %s", sql);
    int64_t v = sqlite3_column_int64(s, 0);
    sqlite3_finalize(s);
    return v;
}

#ifdef ISAM_ENABLE_LMDB
static void lmdb_set_header_meta(lmdb_backend_t *lmdb, sam_hdr_t *hdr) {
    const char *htxt = sam_hdr_str(hdr);
    lmdb_store_put_meta_i64(lmdb, "store_header_hash",
                            header_text_hash(htxt ? htxt : ""));

    kstring_t sq = { 0, 0, NULL };
    int n = sam_hdr_nref(hdr);
    for (int i = 0; i < n; ++i) {
        char lbuf[32];
        snprintf(lbuf, sizeof lbuf, "%" PRId64,
                 (int64_t)sam_hdr_tid2len(hdr, i));
        kputs(sam_hdr_tid2name(hdr, i), &sq);
        kputc('\t', &sq);
        kputs(lbuf, &sq);
        kputc('\n', &sq);
    }
    lmdb_store_put_meta_i64(lmdb, "store_sq_hash",
                            header_text_hash(sq.s ? sq.s : ""));
    free(sq.s);
}

static int index_lmdb(const char *store_path, const char *out_path,
                      const char *name_alias, samFile *sam,
                      sam_hdr_t *hdr, int is_bam,
                      const struct stat *stbuf) {
    struct stat out_st;
    if (stat(out_path, &out_st) == 0)
        die("LMDB output path '%s' already exists; choose a new path", out_path);
    if (errno != ENOENT)
        die("stat '%s': %s", out_path, strerror(errno));

    lmdb_backend_t *lmdb = lmdb_backend_open(out_path, /*readonly*/0, 0);
    lmdb_backend_begin(lmdb);
    lmdb_store_put_meta(lmdb, "backend", "lmdb");
    lmdb_store_put_meta(lmdb, "isam_type", "store");
    lmdb_store_put_meta(lmdb, "isam_version", ISAM_SCHEMA_STORE);
    lmdb_store_put_meta(lmdb, "isamtools_version", ISAM_VERSION);
    lmdb_store_put_meta(lmdb, "store_path", store_path);
    lmdb_store_put_meta(lmdb, "store_format", is_bam ? "bam" : "cram");
    lmdb_store_put_meta_i64(lmdb, "store_size_bytes", (int64_t)stbuf->st_size);
    lmdb_store_put_meta_i64(lmdb, "store_mtime_ns",
                            (int64_t)stbuf->st_mtime * 1000000000LL);
    if (name_alias) lmdb_store_put_meta(lmdb, "name", name_alias);
    char *ts = now_iso8601();
    lmdb_store_put_meta(lmdb, "created_at", ts);
    free(ts);
    lmdb_set_header_meta(lmdb, hdr);

    bam1_t *b = bam_init1();
    int64_t n_records = 0;
    int have_voffset = is_bam;

    while (1) {
        int64_t voffset = -1;
        if (have_voffset)
            voffset = bgzf_tell(sam->fp.bgzf);

        int rc = sam_read1(sam, hdr, b);
        if (rc == -1) break;
        if (rc < -1)
            die("error reading '%s' at record %lld",
                store_path, (long long)n_records);
        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
            continue;

        int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
        lmdb_store_record_t rec = {
            .read_id = n_records + 1,
            .qname = bam_get_qname(b),
            .qname_hash = qname_hash(bam_get_qname(b)),
            .tid = b->core.tid,
            .pos0 = (is_unmapped || b->core.pos < 0) ? -1 : (int64_t)b->core.pos,
            .end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b),
            .flag = b->core.flag,
            .is_unmapped = is_unmapped,
            .bam_voffset = (have_voffset && voffset >= 0) ? voffset : -1,
            .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0
        };
        lmdb_store_put_record(lmdb, &rec);

        ++n_records;
        if (n_records % 100000 == 0) {
            lmdb_backend_commit(lmdb);
            lmdb_backend_begin(lmdb);
            fprintf(stderr, "isamtools index: %lld records...\r",
                    (long long)n_records);
        }
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(sam);

    lmdb_store_put_meta_i64(lmdb, "indexed_name", 1);
    lmdb_store_put_meta_i64(lmdb, "indexed_coord", 1);
    lmdb_backend_commit(lmdb);
    lmdb_backend_close(lmdb);

    fprintf(stderr, "isamtools index: %lld records indexed\n",
            (long long)n_records);
    fprintf(stderr, "isamtools index: wrote LMDB %s\n", out_path);
    return 0;
}
#endif

int cmd_index(int argc, char **argv) {
    const char *name_alias = NULL;
    const char *out_path   = NULL;
    const char *backend    = "sqlite";

    static struct option longopts[] = {
        {"name",   required_argument, 0, 1000},
        {"output", required_argument, 0, 'o'},
        {"backend", required_argument, 0, 1001},
        {"help",   no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: name_alias = optarg; break;
        case 1001: backend    = optarg; break;
        case 'o':  out_path   = optarg; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *store_path = argv[optind];

    char *out_alloc = NULL;
    if (!out_path) {
        out_alloc = isam_index_path_for(store_path);
        out_path = out_alloc;
    }

    samFile *sam = sam_open(store_path, "r");
    if (!sam) die("cannot open '%s'", store_path);
    const htsFormat *fmt = hts_get_format(sam);
    int is_bam = (fmt && fmt->format == bam);
    int is_cram = (fmt && fmt->format == cram);
    if (!is_bam && !is_cram)
        die("unsupported canonical store format (need BAM or CRAM)");

    sam_hdr_t *hdr = sam_hdr_read(sam);
    if (!hdr) die("cannot read header of '%s'", store_path);

    /* metadata */
    struct stat stbuf;
    if (stat(store_path, &stbuf) != 0)
        die("stat '%s': %s", store_path, strerror(errno));

    if (strcmp(backend, "lmdb") == 0) {
#ifdef ISAM_ENABLE_LMDB
        int rc = index_lmdb(store_path, out_path, name_alias, sam, hdr,
                            is_bam, &stbuf);
        free(out_alloc);
        return rc;
#else
        die("LMDB backend requested but this build has no LMDB support");
#endif
    }
    if (strcmp(backend, "sqlite") != 0)
        die("unknown backend '%s' (expected sqlite or lmdb)", backend);

    /* Open output db (create if needed). */
    sqlite3 *db = db_open(out_path, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
    schema_apply_store(db);

    meta_set     (db, "isam_type",          "store");
    meta_set     (db, "isam_version",       ISAM_SCHEMA_STORE);
    meta_set     (db, "isamtools_version",  ISAM_VERSION);
    meta_set     (db, "store_path",         store_path);
    meta_set     (db, "store_format",       is_bam ? "bam" : "cram");
    meta_set_i64 (db, "store_size_bytes",   (int64_t)stbuf.st_size);
    meta_set_i64 (db, "store_mtime_ns",     (int64_t)stbuf.st_mtime * 1000000000LL);
    if (name_alias) meta_set(db, "name", name_alias);
    char *ts = now_iso8601(); meta_set(db, "created_at", ts); free(ts);

    /* hash header text */
    const char *htxt = sam_hdr_str(hdr);
    meta_set_i64(db, "store_header_hash", header_text_hash(htxt ? htxt : ""));
    /* SQ hash: hash only @SQ lines */
    kstring_t sq = { 0, 0, NULL };
    int n = sam_hdr_nref(hdr);
    for (int i = 0; i < n; ++i) {
        char lbuf[32];
        snprintf(lbuf, sizeof lbuf, "%" PRId64,
                 (int64_t)sam_hdr_tid2len(hdr, i));
        kputs(sam_hdr_tid2name(hdr, i), &sq);
        kputc('\t', &sq);
        kputs(lbuf, &sq);
        kputc('\n', &sq);
    }
    meta_set_i64(db, "store_sq_hash", header_text_hash(sq.s ? sq.s : ""));
    free(sq.s);

    /* populate ref_contig */
    db_exec(db, "BEGIN;");
    db_exec(db, "DELETE FROM ref_contig;");
    sqlite3_stmt *ins_ref = db_prep(db,
        "INSERT INTO ref_contig(tid,name,length,md5) VALUES(?,?,?,?);");
    for (int i = 0; i < n; ++i) {
        kstring_t m5 = { 0, 0, NULL };
        const char *md5 = NULL;
        if (sam_hdr_find_tag_pos(hdr, "SQ", i, "M5", &m5) == 0) md5 = m5.s;
        sqlite3_reset(ins_ref);
        sqlite3_bind_int  (ins_ref, 1, i);
        sqlite3_bind_text (ins_ref, 2, sam_hdr_tid2name(hdr,i), -1, SQLITE_STATIC);
        sqlite3_bind_int64(ins_ref, 3, (sqlite3_int64)sam_hdr_tid2len(hdr,i));
        if (md5) sqlite3_bind_text(ins_ref, 4, md5, -1, SQLITE_TRANSIENT);
        else     sqlite3_bind_null(ins_ref, 4);
        db_step_done(ins_ref);
        free(m5.s);
    }
    sqlite3_finalize(ins_ref);
    db_exec(db, "COMMIT;");

    /* Preserve store-level annotations across read_id reassignment. */
    db_exec(db, "DROP TABLE IF EXISTS temp.rebuild_subset_read;");
    db_exec(db, "DROP TABLE IF EXISTS temp.rebuild_read_tag;");
    db_exec(db,
        "CREATE TEMP TABLE rebuild_subset_read AS "
        "SELECT sr.subset_id, cr.qname, cr.qname_hash "
        "FROM subset_read sr "
        "JOIN canonical_read cr ON cr.read_id=sr.read_id;");
    db_exec(db,
        "CREATE TEMP TABLE rebuild_read_tag AS "
        "SELECT cr.qname, cr.qname_hash, rt.tag, rt.type, rt.value, rt.source "
        "FROM read_tag rt "
        "JOIN canonical_read cr ON cr.read_id=rt.read_id;");

    /* clear and re-populate canonical_read */
    db_exec(db, "DELETE FROM read_tag;");
    db_exec(db, "DELETE FROM subset_read;");
    db_exec(db, "DELETE FROM canonical_read;");
    db_exec(db, "BEGIN;");

    sqlite3_stmt *ins = db_prep(db,
        "INSERT INTO canonical_read("
        "  qname, qname_hash, tid, pos0, end0, flag, is_unmapped,"
        "  location_kind, bam_voffset, record_ordinal,"
        "  logical_valid, physical_hint_valid"
        ") VALUES(?,?,?,?,?,?,?,?,?,?,1,?);");

    bam1_t *b = bam_init1();
    int64_t n_records = 0;
    int     have_voffset = is_bam;
    int64_t batch = 0;
    const int64_t BATCH_SIZE = 100000;

    while (1) {
        int64_t voffset = -1;
        if (have_voffset)
            voffset = bgzf_tell(sam->fp.bgzf);

        int rc = sam_read1(sam, hdr, b);
        if (rc == -1) break;
        if (rc < -1)  die("error reading '%s' at record %lld",
                          store_path, (long long)n_records);

        /* Only one canonical record per QNAME: skip secondary/supplementary. */
        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
            continue;

        int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
        int tid  = b->core.tid;
        int64_t pos0 = b->core.pos;
        int64_t end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b);

        sqlite3_reset(ins);
        sqlite3_bind_text  (ins, 1, bam_get_qname(b), -1, SQLITE_TRANSIENT);
        sqlite3_bind_int64 (ins, 2, qname_hash(bam_get_qname(b)));
        if (tid < 0) sqlite3_bind_null(ins, 3); else sqlite3_bind_int(ins, 3, tid);
        if (pos0 < 0 || is_unmapped) sqlite3_bind_null(ins, 4);
        else                          sqlite3_bind_int64(ins, 4, pos0);
        if (end0 < 0) sqlite3_bind_null(ins, 5); else sqlite3_bind_int64(ins, 5, end0);
        sqlite3_bind_int   (ins, 6, b->core.flag);
        sqlite3_bind_int   (ins, 7, is_unmapped);
        sqlite3_bind_text  (ins, 8, is_bam ? "bam_voffset" : "cram_logical",
                            -1, SQLITE_STATIC);
        if (have_voffset && voffset >= 0)
            sqlite3_bind_int64(ins, 9, voffset);
        else
            sqlite3_bind_null (ins, 9);
        sqlite3_bind_int64 (ins, 10, n_records);  /* record_ordinal */
        sqlite3_bind_int   (ins, 11, have_voffset ? 1 : 0);

        db_step_done(ins);

        ++n_records; ++batch;
        if (batch >= BATCH_SIZE) {
            db_exec(db, "COMMIT; BEGIN;");
            batch = 0;
            fprintf(stderr, "isamtools index: %lld records...\r",
                    (long long)n_records);
        }
    }
    db_exec(db, "COMMIT;");
    sqlite3_finalize(ins);
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(sam);

    fprintf(stderr, "isamtools index: %lld records indexed\n",
            (long long)n_records);

    /* build indexes */
    fprintf(stderr, "isamtools index: building indexes...\n");
    schema_build_store_indexes(db);

    int64_t missing_subset = scalar_i64(db,
        "SELECT COUNT(*) FROM rebuild_subset_read o "
        "LEFT JOIN canonical_read cr "
        "  ON cr.qname_hash=o.qname_hash AND cr.qname=o.qname "
        "WHERE cr.read_id IS NULL;");
    int64_t missing_tags = scalar_i64(db,
        "SELECT COUNT(*) FROM rebuild_read_tag o "
        "LEFT JOIN canonical_read cr "
        "  ON cr.qname_hash=o.qname_hash AND cr.qname=o.qname "
        "WHERE cr.read_id IS NULL;");
    if (missing_subset || missing_tags)
        die("store rebuild lost annotated qnames "
            "(%lld subset rows, %lld tag rows missing in new store)",
            (long long)missing_subset, (long long)missing_tags);

    db_exec(db, "BEGIN;");
    db_exec(db,
        "INSERT INTO subset_read(subset_id, read_id) "
        "SELECT o.subset_id, cr.read_id "
        "FROM rebuild_subset_read o "
        "JOIN canonical_read cr "
        "  ON cr.qname_hash=o.qname_hash AND cr.qname=o.qname;");
    db_exec(db,
        "INSERT INTO read_tag(read_id, tag, type, value, source) "
        "SELECT cr.read_id, o.tag, o.type, o.value, o.source "
        "FROM rebuild_read_tag o "
        "JOIN canonical_read cr "
        "  ON cr.qname_hash=o.qname_hash AND cr.qname=o.qname;");
    db_exec(db, "COMMIT;");
    db_exec(db, "DROP TABLE rebuild_subset_read;");
    db_exec(db, "DROP TABLE rebuild_read_tag;");

    meta_set(db, "indexed_name", "1");
    meta_set(db, "indexed_coord", "1");

    db_close(db);
    fprintf(stderr, "isamtools index: wrote %s\n", out_path);
    free(out_alloc);
    return 0;
}
