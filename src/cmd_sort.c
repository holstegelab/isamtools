#include "common.h"
#include "backend_lmdb.h"

#include <getopt.h>
#include <strings.h>
#include <sys/stat.h>
#include <errno.h>

static void usage(void) {
    fputs(
        "Usage: isamtools sort [options] <file.isam|store.ssam|store.bam|store.cram>\n"
        "\n"
        "Build logical indexes for an alignment .isam or canonical store .ssam.\n"
        "\n"
        "Options:\n"
        "  -n             Build name index only.\n"
        "      --coord    Build coordinate-clustered table + indexes (default).\n"
        "      --interval Build LMDB checkpoint interval index.\n"
        "      --all      Build all indexes.\n"
        "  -h, --help     Show this help.\n",
        stderr);
}

/* Resolve user path to an actual LMDB index path. */
static char *resolve_isam_path(const char *p) {
    if (isam_is_view_index_path(p) || isam_is_store_index_path(p))
        return xstrdup(p);
    return isam_store_index_path_for(p);
}

#define LMDB_CHECKPOINT_MIN_SPAN 200
#define LMDB_CHECKPOINT_STEP 16384

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
        int64_t max_ref_span = 0;
        while (lmdb_view_iter_next(it, &rec)) {
            if (!rec.is_unmapped && rec.tid >= 0 && rec.pos0 >= 0) {
                int64_t span = rec.end0 - rec.pos0;
                if (span > max_ref_span) max_ref_span = span;
            }
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
        lmdb_view_put_meta_i64(lmdb, "view_max_ref_span", max_ref_span);
        if (do_coord) {
            lmdb_view_put_meta(lmdb, "indexed_coord", "1");
            lmdb_view_put_meta(lmdb, "storage_layout", "lmdb_record_refs");
        }
        if (do_interval) {
            if (max_ref_span > LMDB_CHECKPOINT_MIN_SPAN) {
                fprintf(stderr, "isamtools sort: building LMDB interval checkpoints...\n");
                lmdb_view_build_checkpoints(lmdb, LMDB_CHECKPOINT_STEP);
                lmdb_view_put_meta(lmdb, "interval_strategy", "checkpoint");
                lmdb_view_put_meta_i64(lmdb, "checkpoint_step", LMDB_CHECKPOINT_STEP);
            } else {
                lmdb_view_put_meta(lmdb, "interval_strategy", "max_span");
            }
            lmdb_view_put_meta(lmdb, "indexed_interval", "1");
        }
    }
    if (do_name) {
        fprintf(stderr, "isamtools sort: marking LMDB name index available...\n");
        lmdb_view_put_meta(lmdb, "indexed_name", "1");
    }
    return 0;
}

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
        die("'%s' is a directory; LMDB index files are file-backed", isam_path);
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
        die("'%s' is not a valid LMDB index (type=%s)",
            isam_path, type ? type : "<missing>");
    }
    free(type);
    lmdb_backend_close(lmdb);
    free(isam_path);
    return rc;
}

int cmd_sort(int argc, char **argv) {
    return cmd_sort_main(argc, argv);
}
