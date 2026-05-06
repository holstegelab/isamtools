#include "common.h"
#include "backend_lmdb.h"
#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>

static void usage(void) {
    fputs(
        "Usage: isamtools update [options] <isam-or-store>\n"
        "\n"
        "Mutate metadata of a store .isam or alignment-view .isam.\n"
        "\n"
        "Options:\n"
        "  --name NAME            Set logical alias.\n"
        "  --store-path PATH      Set canonical store path (alignment view).\n"
        "  --store-index-path P   Set canonical store-index path (alignment view).\n"
        "  --reference-name NAME  Set reference name (store).\n"
        "  --reference-fasta P    Set reference FASTA path (store).\n"
        "  --set KEY=VALUE        Set arbitrary meta key.\n"
        "  --unset KEY            Remove meta key.\n"
        "  -h, --help             Show this help.\n",
        stderr);
}

static char *resolve_isam_path(const char *p) {
    size_t n = strlen(p);
    if (n >= 5 && strcmp(p + n - 5, ".isam") == 0)
        return xstrdup(p);
    return isam_index_path_for(p);
}

int cmd_update(int argc, char **argv) {
    /* collect operations */
    typedef struct { const char *key; const char *value; int unset; } op_t;
    op_t *ops = NULL;
    size_t n_ops = 0, cap = 0;
    #define ADD_OP(k, v, u) do { \
        if (n_ops == cap) { cap = cap ? cap * 2 : 8; ops = xrealloc(ops, cap * sizeof *ops); } \
        ops[n_ops].key = (k); ops[n_ops].value = (v); ops[n_ops].unset = (u); ++n_ops; \
    } while (0)

    static struct option longopts[] = {
        {"name",              required_argument, 0, 1000},
        {"store-path",        required_argument, 0, 1001},
        {"store-index-path",  required_argument, 0, 1002},
        {"reference-name",    required_argument, 0, 1003},
        {"reference-fasta",   required_argument, 0, 1004},
        {"set",               required_argument, 0, 1005},
        {"unset",             required_argument, 0, 1006},
        {"help",              no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: ADD_OP("name",              optarg, 0); break;
        case 1001: ADD_OP("store_path",        optarg, 0); break;
        case 1002: ADD_OP("store_index_path",  optarg, 0); break;
        case 1003: ADD_OP("reference_name",    optarg, 0); break;
        case 1004: ADD_OP("reference_fasta",   optarg, 0); break;
        case 1005: {
            char *eq = strchr(optarg, '=');
            if (!eq) die("--set requires KEY=VALUE");
            *eq = 0;
            ADD_OP(xstrdup(optarg), xstrdup(eq + 1), 0);
            break;
        }
        case 1006: ADD_OP(optarg, NULL, 1); break;
        case 'h':  usage(); free(ops); return 0;
        default:   usage(); free(ops); return 2;
        }
    }
    if (optind >= argc) { usage(); free(ops); return 2; }
    if (n_ops == 0) die("no update operations specified");

    char *isam_path = resolve_isam_path(argv[optind]);

    struct stat st;
    if (stat(isam_path, &st) != 0)
        die("stat '%s': %s", isam_path, strerror(errno));
    if (S_ISDIR(st.st_mode))
        die("'%s' is a directory; LMDB .isam files are file-backed", isam_path);
    if (!isam_file_is_sqlite(isam_path)) {
#ifdef ISAM_ENABLE_LMDB
        lmdb_backend_t *lmdb = lmdb_backend_open(isam_path, /*readonly*/0, 0);
        lmdb_backend_begin(lmdb);
        for (size_t i = 0; i < n_ops; ++i) {
            if (ops[i].unset) lmdb_view_delete_meta(lmdb, ops[i].key);
            else              lmdb_backend_meta_set(lmdb, ops[i].key, ops[i].value);
        }
        lmdb_backend_commit(lmdb);
        lmdb_backend_close(lmdb);
        free(isam_path);
        free(ops);
        return 0;
#else
        die("'%s' is an LMDB .isam, but this build has no LMDB support", isam_path);
#endif
    }

    sqlite3 *db = db_open(isam_path, SQLITE_OPEN_READWRITE);

    db_exec(db, "BEGIN;");
    sqlite3_stmt *del = db_prep(db, "DELETE FROM meta WHERE key=?;");
    for (size_t i = 0; i < n_ops; ++i) {
        if (ops[i].unset) {
            sqlite3_reset(del);
            sqlite3_bind_text(del, 1, ops[i].key, -1, SQLITE_STATIC);
            db_step_done(del);
        } else {
            meta_set(db, ops[i].key, ops[i].value);
        }
    }
    sqlite3_finalize(del);
    db_exec(db, "COMMIT;");

    db_close(db);
    free(isam_path);
    free(ops);
    return 0;
}
