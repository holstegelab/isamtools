#include "common.h"
#include "schema.h"
#include "backend_lmdb.h"

#include <getopt.h>
#include <strings.h>
#include <sys/stat.h>
#include <errno.h>

static void usage(void) {
    fputs(
        "Usage: isamtools sort [options] <file.isam|store.bam|store.cram>\n"
        "\n"
        "Build logical indexes for an alignment .isam (or rebuild store indexes).\n"
        "\n"
        "Options:\n"
        "  -n             Build name index only.\n"
        "      --coord    Build coordinate-clustered table + indexes (default).\n"
        "      --interval Build R*Tree interval index.\n"
        "      --all      Build all indexes.\n"
        "  -h, --help     Show this help.\n",
        stderr);
}

/* Resolve user path to an actual .isam path. */
static char *resolve_isam_path(const char *p) {
    size_t n = strlen(p);
    if (n >= 5 && strcmp(p + n - 5, ".isam") == 0)
        return xstrdup(p);
    return isam_index_path_for(p);
}

static int sort_store(sqlite3 *db) {
    schema_build_store_indexes(db);
    meta_set(db, "indexed_name",  "1");
    meta_set(db, "indexed_coord", "1");
    return 0;
}

static int sort_view(sqlite3 *db,
                     int do_coord, int do_name, int do_interval) {
    if (do_interval) do_coord = 1;
    schema_build_view_clustered(db);

    if (do_coord) {
        fprintf(stderr, "isamtools sort: building coordinate-clustered table...\n");
        db_exec(db, "BEGIN;");
        db_exec(db, "DELETE FROM aln_record;");
        db_exec(db, "DELETE FROM aln_unmapped;");
        db_exec(db, "DELETE FROM aln_id_lookup;");
        db_exec(db,
            "INSERT INTO aln_record("
            "  tid, pos0, record_order, aln_id, store_read_id, end0, flag,"
            "  mapq, n_cigar, cigar, mate_tid, mate_pos0, tlen,"
            "  aux_codec, aux"
            ") "
            "SELECT tid, pos0, record_order, aln_id, store_read_id, end0, flag,"
            "       mapq, n_cigar, cigar, mate_tid, mate_pos0, tlen,"
            "       aux_codec, aux "
            "FROM aln_record_import "
            "WHERE is_unmapped=0 AND tid IS NOT NULL AND pos0 IS NOT NULL "
            "ORDER BY tid, pos0, record_order, aln_id;");
        db_exec(db,
            "INSERT INTO aln_unmapped("
            "  record_order, aln_id, store_read_id, flag, mapq,"
            "  n_cigar, cigar, mate_tid, mate_pos0, tlen, aux_codec, aux"
            ") "
            "SELECT record_order, aln_id, store_read_id, flag, mapq,"
            "       n_cigar, cigar, mate_tid, mate_pos0, tlen, aux_codec, aux "
            "FROM aln_record_import "
            "WHERE is_unmapped=1 OR tid IS NULL OR pos0 IS NULL "
            "ORDER BY record_order, aln_id;");
        db_exec(db,
            "INSERT INTO aln_id_lookup(aln_id, tid, pos0, record_order) "
            "SELECT aln_id, tid, pos0, record_order FROM aln_record;");
        db_exec(db, "COMMIT;");
        meta_set(db, "indexed_coord", "1");
        meta_set(db, "storage_layout", "clustered_without_rowid");
    }
    if (do_name) {
        fprintf(stderr, "isamtools sort: building name index...\n");
        schema_build_view_indexes(db, /*name*/1, /*coord*/0, /*interval*/0);
        meta_set(db, "indexed_name", "1");
    }
    if (do_interval) {
        fprintf(stderr, "isamtools sort: building R*Tree interval index...\n");
        schema_build_view_indexes(db, /*name*/0, /*coord*/0, /*interval*/1);
        db_exec(db, "BEGIN;");
        db_exec(db, "DELETE FROM aln_interval;");
        db_exec(db,
            "INSERT INTO aln_interval(aln_id, tid_min, tid_max, pos_min, pos_max) "
            "SELECT aln_id, tid, tid, pos0, end0 FROM aln_record;");
        db_exec(db, "COMMIT;");
        meta_set(db, "indexed_interval", "1");
    }
    return 0;
}

#ifdef ISAM_ENABLE_LMDB
static int sort_view_lmdb(lmdb_backend_t *lmdb,
                          int do_coord, int do_name, int do_interval) {
    if (do_interval) do_coord = 1;
    if (do_coord || do_interval) {
        fprintf(stderr, "isamtools sort: building LMDB coordinate tables...\n");
        lmdb_backend_begin(lmdb);
        lmdb_view_clear_sorted(lmdb, do_coord, do_interval);
        lmdb_backend_commit(lmdb);

        lmdb_view_iter_t *it = lmdb_view_iter_open(lmdb, LMDB_VIEW_ITER_IMPORT,
                                                   -1, 0, 0);
        lmdb_backend_begin(lmdb);
        lmdb_view_record_t rec;
        int64_t n = 0;
        while (lmdb_view_iter_next(it, &rec)) {
            lmdb_view_index_record(lmdb, &rec, do_coord, do_interval);
            lmdb_view_record_clear(&rec);
            ++n;
            if (n % 50000 == 0) {
                lmdb_backend_commit(lmdb);
                lmdb_backend_begin(lmdb);
                fprintf(stderr, "isamtools sort: %lld records...\r",
                        (long long)n);
            }
        }
        lmdb_backend_commit(lmdb);
        lmdb_view_iter_close(it);
        if (do_coord) {
            lmdb_view_put_meta(lmdb, "indexed_coord", "1");
            lmdb_view_put_meta(lmdb, "storage_layout", "lmdb_record_refs");
        }
        if (do_interval) lmdb_view_put_meta(lmdb, "indexed_interval", "1");
    }
    if (do_name) {
        fprintf(stderr, "isamtools sort: marking LMDB name index available...\n");
        lmdb_view_put_meta(lmdb, "indexed_name", "1");
    }
    return 0;
}
#endif

int cmd_sort_main(int argc, char **argv) {
    int do_name = 0, do_coord = 0, do_interval = 0;
    int explicit = 0;

    static struct option longopts[] = {
        {"coord",    no_argument, 0, 1000},
        {"interval", no_argument, 0, 1001},
        {"all",      no_argument, 0, 1002},
        {"help",     no_argument, 0, 'h'},
        {0,0,0,0}
    };

    /* allow re-entry: reset getopt state */
    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "nh", longopts, NULL)) != -1) {
        switch (c) {
        case 'n':  do_name = 1; explicit = 1; break;
        case 1000: do_coord = 1; explicit = 1; break;
        case 1001: do_interval = 1; explicit = 1; break;
        case 1002: do_name = do_coord = do_interval = 1; explicit = 1; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (!explicit) do_coord = 1;  /* default */

    if (optind >= argc) { usage(); return 2; }
    const char *user_path = argv[optind];
    char *isam_path = resolve_isam_path(user_path);

    struct stat st;
    if (stat(isam_path, &st) != 0)
        die("stat '%s': %s", isam_path, strerror(errno));
    if (S_ISDIR(st.st_mode))
        die("'%s' is a directory; LMDB .isam files are file-backed", isam_path);
    if (!isam_file_is_sqlite(isam_path)) {
#ifdef ISAM_ENABLE_LMDB
        lmdb_backend_t *lmdb = lmdb_backend_open(isam_path, /*readonly*/0, 0);
        char *type = lmdb_backend_meta_get(lmdb, "isam_type");
        int rc;
        if (type && strcmp(type, "store") == 0) {
            lmdb_store_put_meta(lmdb, "indexed_name", "1");
            lmdb_store_put_meta(lmdb, "indexed_coord", "1");
            rc = 0;
        } else if (type && strcmp(type, "alignment") == 0) {
            rc = sort_view_lmdb(lmdb, do_coord, do_name, do_interval);
        } else {
            die("'%s' is not a valid LMDB .isam (type=%s)",
                isam_path, type ? type : "<missing>");
        }
        free(type);
        lmdb_backend_close(lmdb);
        free(isam_path);
        return rc;
#else
        die("'%s' is an LMDB .isam, but this build has no LMDB support", isam_path);
#endif
    }

    sqlite3 *db = db_open(isam_path, SQLITE_OPEN_READWRITE);
    char *type = meta_get(db, "isam_type");
    int rc;
    if (type && strcmp(type, "store") == 0) {
        rc = sort_store(db);
    } else if (type && strcmp(type, "alignment") == 0) {
        rc = sort_view(db, do_coord, do_name, do_interval);
    } else {
        die("'%s' is not a valid .isam (type=%s)",
            isam_path, type ? type : "<missing>");
    }
    free(type);
    db_close(db);
    free(isam_path);
    return rc;
}

int cmd_sort(int argc, char **argv) {
    return cmd_sort_main(argc, argv);
}
