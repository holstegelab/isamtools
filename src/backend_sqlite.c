#include "common.h"
#include "backend.h"

sqlite3 *sqlite_backend_open(const char *path, int flags) {
    sqlite3 *db = NULL;
    int rc = sqlite3_open_v2(path, &db, flags, NULL);
    if (rc != SQLITE_OK) {
        const char *msg = db ? sqlite3_errmsg(db) : sqlite3_errstr(rc);
        die("cannot open sqlite db '%s': %s", path, msg);
    }
    sqlite3_busy_timeout(db, 30000);
    sqlite3_exec(db, "PRAGMA journal_mode=WAL;", NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA synchronous=NORMAL;", NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA temp_store=MEMORY;", NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA foreign_keys=ON;", NULL, NULL, NULL);
    return db;
}

void sqlite_backend_close(sqlite3 *db) {
    if (!db) return;
    int rc = sqlite3_close(db);
    if (rc != SQLITE_OK)
        warnx_("sqlite3_close: %s", sqlite3_errmsg(db));
}

void sqlite_backend_exec(sqlite3 *db, const char *sql) {
    char *err = NULL;
    int rc = sqlite3_exec(db, sql, NULL, NULL, &err);
    if (rc != SQLITE_OK) {
        char *msg = xstrdup(err ? err : sqlite3_errmsg(db));
        sqlite3_free(err);
        die("sqlite exec failed: %s\nSQL: %s", msg, sql);
    }
}

sqlite3_stmt *sqlite_backend_prep(sqlite3 *db, const char *sql) {
    sqlite3_stmt *stmt = NULL;
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK)
        die("sqlite prepare failed: %s\nSQL: %s", sqlite3_errmsg(db), sql);
    return stmt;
}

void sqlite_backend_step_done(sqlite3_stmt *stmt) {
    int rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE)
        die("sqlite step expected DONE, got rc=%d: %s",
            rc, sqlite3_errmsg(sqlite3_db_handle(stmt)));
    sqlite3_reset(stmt);
}

int sqlite_backend_step(sqlite3_stmt *stmt) {
    int rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW && rc != SQLITE_DONE)
        die("sqlite step failed: %s", sqlite3_errmsg(sqlite3_db_handle(stmt)));
    return rc;
}

void backend_meta_set(sqlite3 *db, const char *key, const char *value) {
    sqlite3_stmt *s = sqlite_backend_prep(db,
        "INSERT INTO meta(key,value) VALUES(?,?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value;");
    sqlite3_bind_text(s, 1, key, -1, SQLITE_STATIC);
    sqlite3_bind_text(s, 2, value, -1, SQLITE_TRANSIENT);
    sqlite_backend_step_done(s);
    sqlite3_finalize(s);
}

void backend_meta_set_i64(sqlite3 *db, const char *key, int64_t v) {
    char buf[32];
    snprintf(buf, sizeof buf, "%lld", (long long)v);
    backend_meta_set(db, key, buf);
}

char *backend_meta_get(sqlite3 *db, const char *key) {
    sqlite3_stmt *s = sqlite_backend_prep(db,
        "SELECT value FROM meta WHERE key=?;");
    sqlite3_bind_text(s, 1, key, -1, SQLITE_STATIC);
    char *out = NULL;
    if (sqlite_backend_step(s) == SQLITE_ROW) {
        const unsigned char *t = sqlite3_column_text(s, 0);
        if (t) out = xstrdup((const char *)t);
    }
    sqlite3_finalize(s);
    return out;
}

int backend_meta_get_i64(sqlite3 *db, const char *key, int64_t *out) {
    char *s = backend_meta_get(db, key);
    if (!s) return 0;
    *out = (int64_t)strtoll(s, NULL, 10);
    free(s);
    return 1;
}
