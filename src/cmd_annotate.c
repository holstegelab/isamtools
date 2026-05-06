#include "common.h"
#include "schema.h"

#include <getopt.h>
#include <ctype.h>

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

/* Resolve qname -> integer id (read_id for store, store_read_id for view).
 * Returns 1 on hit, 0 on miss. */
static int lookup_qname(sqlite3_stmt *q, const char *qname, int64_t *out) {
    sqlite3_reset(q);
    sqlite3_bind_int64(q, 1, qname_hash(qname));
    sqlite3_bind_text (q, 2, qname, -1, SQLITE_TRANSIENT);
    int rc = sqlite3_step(q);
    if (rc == SQLITE_ROW) { *out = sqlite3_column_int64(q, 0); return 1; }
    if (rc != SQLITE_DONE)
        die("annotate lookup failed: %s",
            sqlite3_errmsg(sqlite3_db_handle(q)));
    return 0;
}

static void annotate_tags(sqlite3 *db, int is_store,
                          const char *path, const char *cmdline,
                          int allow_missing) {
    FILE *f = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!f) die("cannot open '%s'", path);

    sqlite3_stmt *lookup = db_prep(db, is_store
        ? "SELECT read_id FROM canonical_read "
          "WHERE qname_hash=? AND qname=? LIMIT 1;"
        : "SELECT store_read_id FROM read_dict "
          "WHERE qname_hash=? AND qname=? LIMIT 1;");
    sqlite3_stmt *ins = db_prep(db, is_store
        ? "INSERT INTO read_tag(read_id, tag, type, value, source) "
          "VALUES(?,?,?,?,?) "
          "ON CONFLICT(read_id, tag) DO UPDATE SET "
          "  type=excluded.type, value=excluded.value, source=excluded.source;"
        : "INSERT INTO read_tag(store_read_id, tag, type, value, source) "
          "VALUES(?,?,?,?,?) "
          "ON CONFLICT(store_read_id, tag) DO UPDATE SET "
          "  type=excluded.type, value=excluded.value, source=excluded.source;");

    db_exec(db, "BEGIN;");

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
        if (!lookup_qname(lookup, qname, &id)) {
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

            sqlite3_reset(ins);
            sqlite3_bind_int64(ins, 1, id);
            sqlite3_bind_text (ins, 2, tag,  -1, SQLITE_TRANSIENT);
            sqlite3_bind_text (ins, 3, type, -1, SQLITE_TRANSIENT);
            sqlite3_bind_blob (ins, 4, val, (int)strlen(val), SQLITE_TRANSIENT);
            sqlite3_bind_text (ins, 5, path, -1, SQLITE_STATIC);
            db_step_done(ins);
            ++n_tags;
            cur = next ? next + 1 : NULL;
        }
        ++n_rows;
    }
    db_exec(db, "COMMIT;");
    sqlite3_finalize(ins);
    sqlite3_finalize(lookup);
    free(line);
    if (f != stdin) fclose(f);

    sqlite3_stmt *am = db_prep(db,
        "INSERT INTO annotation_meta(kind, source, command_line, created_at) "
        "VALUES('tag', ?, ?, ?);");
    sqlite3_bind_text(am, 1, path,    -1, SQLITE_STATIC);
    sqlite3_bind_text(am, 2, cmdline, -1, SQLITE_STATIC);
    char *ts = now_iso8601();
    sqlite3_bind_text(am, 3, ts, -1, SQLITE_TRANSIENT);
    db_step_done(am);
    sqlite3_finalize(am);
    free(ts);

    fprintf(stderr,
            "isamtools annotate: %lld rows, %lld tags, %lld missing\n",
            (long long)n_rows, (long long)n_tags, (long long)n_missing);
}

static void annotate_subset(sqlite3 *db, int is_store,
                            const char *subset_name, const char *desc,
                            const char *path, const char *cmdline,
                            int allow_missing) {
    sqlite3_stmt *up = db_prep(db,
        "INSERT INTO subset(name, description, created_at) VALUES(?,?,?) "
        "ON CONFLICT(name) DO UPDATE SET "
        "  description=COALESCE(excluded.description, subset.description) "
        "RETURNING subset_id;");
    sqlite3_bind_text(up, 1, subset_name, -1, SQLITE_STATIC);
    if (desc) sqlite3_bind_text(up, 2, desc, -1, SQLITE_STATIC);
    else      sqlite3_bind_null(up, 2);
    char *ts = now_iso8601();
    sqlite3_bind_text(up, 3, ts, -1, SQLITE_TRANSIENT);
    int rc = sqlite3_step(up);
    if (rc != SQLITE_ROW) die("subset upsert failed: %s", sqlite3_errmsg(db));
    int64_t subset_id = sqlite3_column_int64(up, 0);
    sqlite3_finalize(up);
    free(ts);

    sqlite3_stmt *lookup = db_prep(db, is_store
        ? "SELECT read_id FROM canonical_read "
          "WHERE qname_hash=? AND qname=? LIMIT 1;"
        : "SELECT store_read_id FROM read_dict "
          "WHERE qname_hash=? AND qname=? LIMIT 1;");
    sqlite3_stmt *ins = db_prep(db, is_store
        ? "INSERT OR IGNORE INTO subset_read(subset_id, read_id) VALUES(?,?);"
        : "INSERT OR IGNORE INTO subset_read(subset_id, store_read_id) "
          "VALUES(?,?);");

    FILE *f = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!f) die("cannot open '%s'", path);

    db_exec(db, "BEGIN;");
    char *line = NULL; size_t cap = 0;
    ssize_t n;
    int64_t n_added = 0, n_missing = 0;
    while ((n = getline(&line, &cap, f)) >= 0) {
        rstrip(line);
        if (!*line || *line == '#') continue;

        int64_t id;
        if (!lookup_qname(lookup, line, &id)) {
            if (!allow_missing)
                die("annotate: qname '%s' not present in target "
                    "(use --allow-missing to skip)", line);
            ++n_missing;
            continue;
        }

        sqlite3_reset(ins);
        sqlite3_bind_int64(ins, 1, subset_id);
        sqlite3_bind_int64(ins, 2, id);
        db_step_done(ins);
        ++n_added;
    }
    db_exec(db, "COMMIT;");
    sqlite3_finalize(ins);
    sqlite3_finalize(lookup);
    free(line);
    if (f != stdin) fclose(f);

    sqlite3_stmt *am = db_prep(db,
        "INSERT INTO annotation_meta(kind, source, command_line, created_at) "
        "VALUES('subset', ?, ?, ?);");
    sqlite3_bind_text(am, 1, path,    -1, SQLITE_STATIC);
    sqlite3_bind_text(am, 2, cmdline, -1, SQLITE_STATIC);
    char *ts2 = now_iso8601();
    sqlite3_bind_text(am, 3, ts2, -1, SQLITE_TRANSIENT);
    db_step_done(am);
    sqlite3_finalize(am);
    free(ts2);

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
    sqlite3 *db = db_open(isam_path, SQLITE_OPEN_READWRITE);
    char *type = meta_get(db, "isam_type");
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

    if (subset) annotate_subset(db, is_store, subset, desc, from, cmdbuf,
                                allow_missing);
    else        annotate_tags  (db, is_store, from, cmdbuf, allow_missing);

    db_close(db);
    free(isam_path);
    return 0;
}
