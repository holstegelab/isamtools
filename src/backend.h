#ifndef ISAM_BACKEND_H
#define ISAM_BACKEND_H

#include <stdint.h>
#include <sqlite3.h>

typedef sqlite3 isam_sqlite_db;
typedef sqlite3_stmt isam_sqlite_stmt;

sqlite3      *sqlite_backend_open(const char *path, int flags);
void          sqlite_backend_close(sqlite3 *db);
void          sqlite_backend_exec(sqlite3 *db, const char *sql);
sqlite3_stmt *sqlite_backend_prep(sqlite3 *db, const char *sql);
void          sqlite_backend_step_done(sqlite3_stmt *stmt);
int           sqlite_backend_step(sqlite3_stmt *stmt);

void  backend_meta_set(sqlite3 *db, const char *key, const char *value);
void  backend_meta_set_i64(sqlite3 *db, const char *key, int64_t v);
char *backend_meta_get(sqlite3 *db, const char *key);
int   backend_meta_get_i64(sqlite3 *db, const char *key, int64_t *out);

/* Compatibility aliases while command code is migrated to semantic backends. */
#define db_open       sqlite_backend_open
#define db_close      sqlite_backend_close
#define db_exec       sqlite_backend_exec
#define db_prep       sqlite_backend_prep
#define db_step_done  sqlite_backend_step_done
#define db_step       sqlite_backend_step

#define meta_set      backend_meta_set
#define meta_set_i64  backend_meta_set_i64
#define meta_get      backend_meta_get
#define meta_get_i64  backend_meta_get_i64

#endif
