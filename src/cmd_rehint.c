#include "common.h"
#include "backend_lmdb.h"

#include <getopt.h>

static void usage(void) {
    fputs(
        "Usage: isamtools rehint [options] <view.isam>\n"
        "\n"
        "Re-resolve every alignment-view store_read_id through the current\n"
        "canonical store index using read_dict qnames.\n"
        "\n"
        "Options:\n"
        "  --store PATH         Update/use canonical store path metadata.\n"
        "  --store-index PATH   Use canonical store .ssam (default from view).\n"
        "  -h, --help           Show this help.\n",
        stderr);
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
        case 'h': usage(); return 0;
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *view_path = argv[optind];

    lmdb_backend_t *view = lmdb_backend_open(view_path, /*readonly*/0, 0);
    char *view_type = lmdb_backend_meta_get(view, "isam_type");
    if (!view_type || strcmp(view_type, "alignment") != 0)
        die("'%s' is not an alignment-view .isam (type=%s)",
            view_path, view_type ? view_type : "<missing>");
    free(view_type);

    char *store_path_meta = lmdb_backend_meta_get(view, "store_path");
    char *store_idx_meta = lmdb_backend_meta_get(view, "store_index_path");
    char *resolved_path = isam_path_relative_to_index(view_path, store_path_meta);
    char *resolved_idx = isam_path_relative_to_index(view_path, store_idx_meta);
    const char *store_path = store_path_override ? store_path_override
                                                 : resolved_path;
    const char *store_idx = store_idx_override ? store_idx_override
                                               : resolved_idx;
    char *store_idx_alloc = NULL;
    if (!store_idx && store_path) {
        store_idx_alloc = isam_store_index_path_for(store_path);
        store_idx = store_idx_alloc;
    }
    if (!store_idx)
        die("rehint: no store index path available; pass --store-index");

    lmdb_backend_t *store = lmdb_backend_open(store_idx, /*readonly*/1, 0);
    char *store_type = lmdb_backend_meta_get(store, "isam_type");
    if (!store_type || strcmp(store_type, "store") != 0)
        die("'%s' is not a store index .ssam (type=%s)",
            store_idx, store_type ? store_type : "<missing>");
    free(store_type);

    int64_t n_reads = 0, n_changed = 0;
    /* Preserve portable metadata unless the caller replaces it. Explicit
     * overrides must persist their cwd-relative meaning as absolute paths. */
    char *persisted_path = store_path_override ? isam_path_absolute(store_path_override)
                                               : xstrdup(store_path_meta);
    char *persisted_idx = store_idx_override ? isam_path_absolute(store_idx_override)
                            : store_idx_meta ? xstrdup(store_idx_meta)
                                             : isam_path_absolute(store_idx);
    lmdb_view_rehint(view, store, persisted_path, persisted_idx, &n_reads, &n_changed);

    lmdb_backend_close(store);
    lmdb_backend_close(view);
    free(store_idx_alloc);
    free(persisted_path);
    free(persisted_idx);
    free(resolved_path);
    free(resolved_idx);
    free(store_path_meta);
    free(store_idx_meta);

    fprintf(stderr,
            "isamtools rehint: %lld reads resolved, %lld ids changed\n",
            (long long)n_reads, (long long)n_changed);
    return 0;
}
