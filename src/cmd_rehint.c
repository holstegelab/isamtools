#include "common.h"

#include <getopt.h>

static void usage(void) {
    fputs(
        "Usage: isamtools rehint [options] <view.isam>\n"
        "\n"
        "Re-resolve every alignment-view store_read_id through the current\n"
        "canonical store index using read_dict.qname.\n"
        "\n"
        "Options:\n"
        "  --store PATH         Update/use canonical store path metadata.\n"
        "  --store-index PATH   Use canonical store .isam (default from view).\n"
        "  -h, --help           Show this help.\n",
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

static int table_exists(sqlite3 *db, const char *name) {
    sqlite3_stmt *s = db_prep(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? "
        "UNION ALL "
        "SELECT 1 FROM sqlite_temp_master WHERE type='table' AND name=? "
        "LIMIT 1;");
    sqlite3_bind_text(s, 1, name, -1, SQLITE_STATIC);
    sqlite3_bind_text(s, 2, name, -1, SQLITE_STATIC);
    int exists = sqlite3_step(s) == SQLITE_ROW;
    sqlite3_finalize(s);
    return exists;
}

static void require_all_mapped(sqlite3 *db, const char *table) {
    char sql[512];
    snprintf(sql, sizeof sql,
        "SELECT COUNT(*) FROM %s t "
        "LEFT JOIN rehint_map m ON m.old_store_read_id=t.store_read_id "
        "WHERE m.new_store_read_id IS NULL;",
        table);
    int64_t missing = scalar_i64(db, sql);
    if (missing)
        die("rehint: table %s contains %lld store_read_id values not present "
            "in read_dict", table, (long long)missing);
}

static void rewrite_aln_record_import(sqlite3 *db) {
    require_all_mapped(db, "aln_record_import");
    db_exec(db,
        "CREATE TEMP TABLE rehint_aln_record_import AS "
        "SELECT t.aln_id, m.new_store_read_id AS store_read_id,"
        "       t.flag, t.tid, t.pos0, t.end0, t.mapq, t.n_cigar, t.cigar,"
        "       t.mate_tid, t.mate_pos0, t.tlen, t.aux_codec, t.aux,"
        "       t.record_order, t.is_unmapped "
        "FROM aln_record_import t "
        "JOIN rehint_map m ON m.old_store_read_id=t.store_read_id;");
    db_exec(db, "DELETE FROM aln_record_import;");
    db_exec(db,
        "INSERT INTO aln_record_import("
        "  aln_id, store_read_id, flag, tid, pos0, end0, mapq, n_cigar,"
        "  cigar, mate_tid, mate_pos0, tlen, aux_codec, aux, record_order,"
        "  is_unmapped"
        ") SELECT aln_id, store_read_id, flag, tid, pos0, end0, mapq,"
        "         n_cigar, cigar, mate_tid, mate_pos0, tlen, aux_codec, aux,"
        "         record_order, is_unmapped "
        "FROM rehint_aln_record_import;");
    db_exec(db, "DROP TABLE rehint_aln_record_import;");
}

static void rewrite_aln_record(sqlite3 *db) {
    if (!table_exists(db, "aln_record")) return;
    require_all_mapped(db, "aln_record");
    db_exec(db,
        "CREATE TEMP TABLE rehint_aln_record AS "
        "SELECT t.tid, t.pos0, t.record_order, t.aln_id,"
        "       m.new_store_read_id AS store_read_id, t.end0, t.flag,"
        "       t.mapq, t.n_cigar, t.cigar, t.mate_tid, t.mate_pos0,"
        "       t.tlen, t.aux_codec, t.aux "
        "FROM aln_record t "
        "JOIN rehint_map m ON m.old_store_read_id=t.store_read_id;");
    db_exec(db, "DELETE FROM aln_record;");
    db_exec(db,
        "INSERT INTO aln_record("
        "  tid, pos0, record_order, aln_id, store_read_id, end0, flag,"
        "  mapq, n_cigar, cigar, mate_tid, mate_pos0, tlen, aux_codec, aux"
        ") SELECT tid, pos0, record_order, aln_id, store_read_id, end0,"
        "         flag, mapq, n_cigar, cigar, mate_tid, mate_pos0, tlen,"
        "         aux_codec, aux "
        "FROM rehint_aln_record;");
    db_exec(db, "DROP TABLE rehint_aln_record;");
}

static void rewrite_aln_unmapped(sqlite3 *db) {
    if (!table_exists(db, "aln_unmapped")) return;
    require_all_mapped(db, "aln_unmapped");
    db_exec(db,
        "CREATE TEMP TABLE rehint_aln_unmapped AS "
        "SELECT t.record_order, t.aln_id, m.new_store_read_id AS store_read_id,"
        "       t.flag, t.mapq, t.n_cigar, t.cigar, t.mate_tid, t.mate_pos0,"
        "       t.tlen, t.aux_codec, t.aux "
        "FROM aln_unmapped t "
        "JOIN rehint_map m ON m.old_store_read_id=t.store_read_id;");
    db_exec(db, "DELETE FROM aln_unmapped;");
    db_exec(db,
        "INSERT INTO aln_unmapped("
        "  record_order, aln_id, store_read_id, flag, mapq, n_cigar, cigar,"
        "  mate_tid, mate_pos0, tlen, aux_codec, aux"
        ") SELECT record_order, aln_id, store_read_id, flag, mapq, n_cigar,"
        "         cigar, mate_tid, mate_pos0, tlen, aux_codec, aux "
        "FROM rehint_aln_unmapped;");
    db_exec(db, "DROP TABLE rehint_aln_unmapped;");
}

static void rewrite_subset_read(sqlite3 *db) {
    require_all_mapped(db, "subset_read");
    db_exec(db,
        "CREATE TEMP TABLE rehint_subset_read AS "
        "SELECT t.subset_id, m.new_store_read_id AS store_read_id "
        "FROM subset_read t "
        "JOIN rehint_map m ON m.old_store_read_id=t.store_read_id;");
    db_exec(db, "DELETE FROM subset_read;");
    db_exec(db,
        "INSERT INTO subset_read(subset_id, store_read_id) "
        "SELECT subset_id, store_read_id FROM rehint_subset_read;");
    db_exec(db, "DROP TABLE rehint_subset_read;");
}

static void rewrite_read_tag(sqlite3 *db) {
    require_all_mapped(db, "read_tag");
    db_exec(db,
        "CREATE TEMP TABLE rehint_read_tag AS "
        "SELECT m.new_store_read_id AS store_read_id, t.tag, t.type,"
        "       t.value, t.source "
        "FROM read_tag t "
        "JOIN rehint_map m ON m.old_store_read_id=t.store_read_id;");
    db_exec(db, "DELETE FROM read_tag;");
    db_exec(db,
        "INSERT INTO read_tag(store_read_id, tag, type, value, source) "
        "SELECT store_read_id, tag, type, value, source FROM rehint_read_tag;");
    db_exec(db, "DROP TABLE rehint_read_tag;");
}

static void rewrite_read_dict(sqlite3 *db) {
    db_exec(db,
        "CREATE TEMP TABLE rehint_read_dict AS "
        "SELECT m.new_store_read_id AS store_read_id, d.qname, d.qname_hash "
        "FROM read_dict d "
        "JOIN rehint_map m ON m.old_store_read_id=d.store_read_id;");
    db_exec(db, "DELETE FROM read_dict;");
    db_exec(db,
        "INSERT INTO read_dict(store_read_id, qname, qname_hash) "
        "SELECT store_read_id, qname, qname_hash FROM rehint_read_dict;");
    db_exec(db, "DROP TABLE rehint_read_dict;");
}

int cmd_rehint(int argc, char **argv) {
    const char *store_path_override = NULL;
    const char *store_idx_override = NULL;

    static struct option longopts[] = {
        {"store",       required_argument, 0, 1000},
        {"store-index", required_argument, 0, 1001},
        {"help",        no_argument,       0, 'h'},
        {0,0,0,0}
    };

    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: store_path_override = optarg; break;
        case 1001: store_idx_override = optarg; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *view_path = argv[optind];

    sqlite3 *view_db = db_open(view_path, SQLITE_OPEN_READWRITE);
    char *view_type = meta_get(view_db, "isam_type");
    if (!view_type || strcmp(view_type, "alignment") != 0)
        die("'%s' is not an alignment-view .isam (type=%s)",
            view_path, view_type ? view_type : "<missing>");
    free(view_type);

    char *store_path_meta = meta_get(view_db, "store_path");
    char *store_idx_meta = meta_get(view_db, "store_index_path");
    const char *store_path = store_path_override ? store_path_override
                                                 : store_path_meta;
    const char *store_idx = store_idx_override ? store_idx_override
                                               : store_idx_meta;
    char *store_idx_alloc = NULL;
    if (!store_idx && store_path) {
        store_idx_alloc = isam_index_path_for(store_path);
        store_idx = store_idx_alloc;
    }
    if (!store_idx)
        die("rehint: no store index path available; pass --store-index");

    sqlite3 *store_db = db_open(store_idx, SQLITE_OPEN_READONLY);
    char *store_type = meta_get(store_db, "isam_type");
    if (!store_type || strcmp(store_type, "store") != 0)
        die("'%s' is not a store index .isam (type=%s)",
            store_idx, store_type ? store_type : "<missing>");
    free(store_type);

    sqlite3_stmt *reads = db_prep(view_db,
        "SELECT store_read_id, qname, qname_hash FROM read_dict "
        "ORDER BY store_read_id;");
    sqlite3_stmt *lookup = db_prep(store_db,
        "SELECT read_id FROM canonical_read "
        "WHERE qname_hash=? AND qname=? LIMIT 1;");
    sqlite3_stmt *ins_map;

    db_exec(view_db, "DROP TABLE IF EXISTS temp.rehint_map;");
    db_exec(view_db,
        "CREATE TEMP TABLE rehint_map("
        "  old_store_read_id INTEGER PRIMARY KEY,"
        "  new_store_read_id INTEGER NOT NULL UNIQUE,"
        "  qname TEXT NOT NULL"
        ");");
    ins_map = db_prep(view_db,
        "INSERT INTO rehint_map(old_store_read_id, new_store_read_id, qname) "
        "VALUES(?,?,?);");

    int64_t n_reads = 0, n_changed = 0;
    while (db_step(reads) == SQLITE_ROW) {
        int64_t old_id = sqlite3_column_int64(reads, 0);
        const char *qname = (const char *)sqlite3_column_text(reads, 1);
        int64_t qh = sqlite3_column_int64(reads, 2);

        sqlite3_reset(lookup);
        sqlite3_bind_int64(lookup, 1, qh);
        sqlite3_bind_text (lookup, 2, qname, -1, SQLITE_TRANSIENT);
        int rc = sqlite3_step(lookup);
        if (rc == SQLITE_DONE)
            die("rehint: qname '%s' is absent from store index '%s'",
                qname, store_idx);
        if (rc != SQLITE_ROW)
            die("rehint lookup failed: %s", sqlite3_errmsg(store_db));
        int64_t new_id = sqlite3_column_int64(lookup, 0);

        sqlite3_reset(ins_map);
        sqlite3_bind_int64(ins_map, 1, old_id);
        sqlite3_bind_int64(ins_map, 2, new_id);
        sqlite3_bind_text (ins_map, 3, qname, -1, SQLITE_TRANSIENT);
        db_step_done(ins_map);

        if (old_id != new_id) ++n_changed;
        ++n_reads;
    }
    sqlite3_finalize(ins_map);
    sqlite3_finalize(lookup);
    sqlite3_finalize(reads);

    int64_t mapped = scalar_i64(view_db, "SELECT COUNT(*) FROM rehint_map;");
    if (mapped != n_reads)
        die("rehint: internal mapping count mismatch");

    db_exec(view_db, "PRAGMA foreign_keys=OFF;");
    db_exec(view_db, "BEGIN;");
    rewrite_aln_record_import(view_db);
    rewrite_aln_record(view_db);
    rewrite_aln_unmapped(view_db);
    rewrite_subset_read(view_db);
    rewrite_read_tag(view_db);
    rewrite_read_dict(view_db);

    if (store_path) meta_set(view_db, "store_path", store_path);
    meta_set(view_db, "store_index_path", store_idx);
    char *shh = meta_get(store_db, "store_header_hash");
    if (shh) { meta_set(view_db, "store_header_hash", shh); free(shh); }
    char *ssh = meta_get(store_db, "store_sq_hash");
    if (ssh) { meta_set(view_db, "store_sq_hash", ssh); free(ssh); }
    char *ts = now_iso8601();
    meta_set(view_db, "rehinted_at", ts);
    free(ts);
    meta_set_i64(view_db, "rehinted_read_count", n_reads);
    db_exec(view_db, "COMMIT;");
    db_exec(view_db, "PRAGMA foreign_keys=ON;");

    db_close(store_db);
    db_close(view_db);
    free(store_idx_alloc);
    free(store_path_meta);
    free(store_idx_meta);

    fprintf(stderr,
            "isamtools rehint: %lld reads resolved, %lld ids changed\n",
            (long long)n_reads, (long long)n_changed);
    return 0;
}
