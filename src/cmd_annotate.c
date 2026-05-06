#include "common.h"
#include "backend_lmdb.h"

#include <getopt.h>

static void usage(void) {
    fputs(
        "Usage: isamtools annotate [options] <isam-or-store>\n"
        "\n"
        "Annotate a store .isam or alignment-view .isam.\n"
        "\n"
        "Tag mode (default):\n"
        "  --from FILE     SAM-tag-like TSV: qname<TAB>TAG:TYPE:VALUE...\n"
        "\n"
        "Subset mode:\n"
        "  --subset NAME   Subset name to add reads to.\n"
        "  --from FILE     One QNAME per line.\n"
        "  --description S Optional subset description.\n"
        "\n"
        "Options:\n"
        "  --allow-missing Skip qnames not present in the target .isam.\n"
        "  -h, --help      Show this help.\n",
        stderr);
}

static char *resolve_isam_path(const char *p) {
    size_t n = strlen(p);
    if (n >= 5 && strcmp(p + n - 5, ".isam") == 0)
        return xstrdup(p);
    return isam_index_path_for(p);
}

static void rstrip(char *s) {
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == '\n' || s[n-1] == '\r' || s[n-1] == ' '
                                    || s[n-1] == '\t'))
        s[--n] = 0;
}

static int lookup_qname(lmdb_backend_t *lmdb,
                        int is_store,
                        const char *qname,
                        int64_t *out) {
    if (is_store) {
        lmdb_store_record_t rec;
        if (!lmdb_store_lookup_qname(lmdb, qname_hash(qname), qname, &rec))
            return 0;
        *out = rec.read_id;
        lmdb_store_record_clear(&rec);
        return 1;
    }
    return lmdb_view_lookup_qname(lmdb, qname_hash(qname), qname, out);
}

static void annotate_tags(lmdb_backend_t *lmdb, int is_store,
                          const char *path, const char *cmdline,
                          int allow_missing) {
    FILE *f = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!f) die("cannot open '%s'", path);

    char *line = NULL; size_t cap = 0;
    ssize_t n;
    int64_t n_rows = 0, n_tags = 0, n_missing = 0;
    while ((n = getline(&line, &cap, f)) >= 0) {
        rstrip(line);
        if (!*line || *line == '#') continue;

        char *qname = line;
        char *tab = strchr(line, '\t');
        if (!tab) die("annotate: malformed line (no TAB): '%s'", line);
        *tab = 0;

        int64_t id;
        if (!lookup_qname(lmdb, is_store, qname, &id)) {
            if (!allow_missing)
                die("annotate: qname '%s' not present in target "
                    "(use --allow-missing to skip)", qname);
            ++n_missing;
            continue;
        }

        char *cur = tab + 1;
        while (cur && *cur) {
            char *next = strchr(cur, '\t');
            if (next) *next = 0;
            if (strlen(cur) < 5 || cur[2] != ':' || cur[4] != ':')
                die("annotate: malformed tag '%s' (need TAG:TYPE:VALUE)", cur);
            char tag[3]  = { cur[0], cur[1], 0 };
            char type[2] = { cur[3], 0 };
            const char *val = cur + 5;

            lmdb_read_tag_put(lmdb, id, tag, type, val, strlen(val), path);
            ++n_tags;
            cur = next ? next + 1 : NULL;
        }
        ++n_rows;
    }
    lmdb_annotation_meta_put(lmdb, "tag", path, cmdline);
    free(line);
    if (f != stdin) fclose(f);

    fprintf(stderr,
            "isamtools annotate: %lld rows, %lld tags, %lld missing\n",
            (long long)n_rows, (long long)n_tags, (long long)n_missing);
}

static void annotate_subset(lmdb_backend_t *lmdb, int is_store,
                            const char *subset_name, const char *desc,
                            const char *path, const char *cmdline,
                            int allow_missing) {
    FILE *f = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!f) die("cannot open '%s'", path);

    int64_t subset_id = lmdb_subset_upsert(lmdb, subset_name, desc);
    char *line = NULL; size_t cap = 0;
    ssize_t n;
    int64_t n_added = 0, n_missing = 0;
    while ((n = getline(&line, &cap, f)) >= 0) {
        rstrip(line);
        if (!*line || *line == '#') continue;

        int64_t id;
        if (!lookup_qname(lmdb, is_store, line, &id)) {
            if (!allow_missing)
                die("annotate: qname '%s' not present in target "
                    "(use --allow-missing to skip)", line);
            ++n_missing;
            continue;
        }

        lmdb_subset_add_read(lmdb, subset_id, id);
        ++n_added;
    }
    lmdb_annotation_meta_put(lmdb, "subset", path, cmdline);
    free(line);
    if (f != stdin) fclose(f);

    fprintf(stderr, "isamtools annotate: subset '%s' += %lld qnames, %lld missing\n",
            subset_name, (long long)n_added, (long long)n_missing);
}

int cmd_annotate(int argc, char **argv) {
    const char *from   = NULL;
    const char *subset = NULL;
    const char *desc   = NULL;
    int allow_missing  = 0;

    static struct option longopts[] = {
        {"from",          required_argument, 0, 1000},
        {"subset",        required_argument, 0, 1001},
        {"description",   required_argument, 0, 1002},
        {"allow-missing", no_argument,       0, 1003},
        {"help",          no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: from   = optarg; break;
        case 1001: subset = optarg; break;
        case 1002: desc   = optarg; break;
        case 1003: allow_missing = 1; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    if (!from) die("--from is required");

    char *isam_path = resolve_isam_path(argv[optind]);
    lmdb_backend_t *lmdb = lmdb_backend_open(isam_path, /*readonly*/0, 0);
    char *type = lmdb_backend_meta_get(lmdb, "isam_type");
    if (!type) die("'%s' has no isam_type meta", isam_path);
    int is_store = (strcmp(type, "store") == 0);
    if (!is_store && strcmp(type, "alignment") != 0)
        die("'%s' is not a recognized .isam (type=%s)", isam_path, type);
    free(type);

    char cmdbuf[1024]; cmdbuf[0] = 0;
    for (int i = 0; i < argc; ++i) {
        size_t cur = strlen(cmdbuf);
        if (cur >= sizeof cmdbuf - 2) break;
        snprintf(cmdbuf + cur, sizeof cmdbuf - cur,
                 "%s%s", i ? " " : "", argv[i]);
    }

    if (subset) annotate_subset(lmdb, is_store, subset, desc, from, cmdbuf,
                                allow_missing);
    else        annotate_tags  (lmdb, is_store, from, cmdbuf, allow_missing);

    lmdb_backend_close(lmdb);
    free(isam_path);
    return 0;
}
