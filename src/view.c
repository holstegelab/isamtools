#include "view.h"
#include "schema.h"

#include <htslib/kstring.h>
#include <sys/stat.h>
#include <errno.h>

view_t *view_open_read(const char *path) {
    view_t *v = xcalloc(1, sizeof *v);
    v->path = xstrdup(path);

    struct stat st;
    if (stat(path, &st) != 0)
        die("stat view '%s' failed: %s", path, strerror(errno));
    if (S_ISDIR(st.st_mode))
        die("view '%s' is a directory; LMDB views are file-backed .isam files",
            path);
    if (!isam_file_is_sqlite(path)) {
#ifdef ISAM_ENABLE_LMDB
        v->backend_kind = VIEW_BACKEND_LMDB;
        v->lmdb = lmdb_backend_open(path, /*readonly*/1, 0);
#else
        die("'%s' is an LMDB alignment view, but this build has no LMDB support",
            path);
#endif
    } else {
        v->backend_kind = VIEW_BACKEND_SQLITE;
        v->db = db_open(path, SQLITE_OPEN_READONLY);
    }

    char *type = view_meta_get(v, "isam_type");
    if (!type || strcmp(type, "alignment") != 0)
        die("'%s' is not an alignment-view .isam (type=%s)",
            path, type ? type : "<missing>");
    free(type);

    v->hdr = view_load_header(v);
    v->n_targets = v->hdr ? sam_hdr_nref(v->hdr) : 0;
    return v;
}

view_t *view_open_write(const char *path) {
    return view_open_write_backend(path, VIEW_BACKEND_SQLITE);
}

view_t *view_open_write_backend(const char *path, view_backend_kind_t backend) {
    view_t *v = xcalloc(1, sizeof *v);
    v->path = xstrdup(path);
    v->backend_kind = backend;
    if (backend == VIEW_BACKEND_SQLITE) {
        v->db = db_open(path, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
        schema_apply_view(v->db);
    } else {
#ifdef ISAM_ENABLE_LMDB
        struct stat st;
        if (stat(path, &st) == 0)
            die("LMDB output path '%s' already exists; choose a new path", path);
        if (errno != ENOENT)
            die("stat '%s': %s", path, strerror(errno));
        v->lmdb = lmdb_backend_open(path, /*readonly*/0, 0);
#else
        die("this build has no LMDB support");
#endif
    }
    return v;
}

void view_close(view_t *v) {
    if (!v) return;
    if (v->store) store_close(v->store);
    if (v->hdr)   sam_hdr_destroy(v->hdr);
#ifdef ISAM_ENABLE_LMDB
    if (v->lmdb)  lmdb_backend_close(v->lmdb);
#endif
    if (v->db)    db_close(v->db);
    free(v->path);
    free(v);
}

void view_set_header(view_t *v, const sam_hdr_t *hdr) {
    /* persist header text */
    const char *txt = sam_hdr_str((sam_hdr_t *)hdr);
    size_t      tln = sam_hdr_length((sam_hdr_t *)hdr);
    if (!txt) die("sam_hdr_str returned NULL");

    if (v->backend_kind == VIEW_BACKEND_LMDB) {
#ifdef ISAM_ENABLE_LMDB
        lmdb_view_put_meta(v->lmdb, "alignment_header_text", txt);
        lmdb_view_put_meta_i64(v->lmdb, "alignment_header_hash",
                               header_text_hash(txt));
        return;
#else
        die("this build has no LMDB support");
#endif
    }

    sqlite3_stmt *s = db_prep(v->db,
        "INSERT INTO meta(key,value) VALUES('alignment_header_text',?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value;");
    sqlite3_bind_text(s, 1, txt, (int)tln, SQLITE_TRANSIENT);
    db_step_done(s);
    sqlite3_finalize(s);

    meta_set_i64(v->db, "alignment_header_hash", header_text_hash(txt));

    /* populate ref_contig */
    db_exec(v->db, "DELETE FROM ref_contig;");
    sqlite3_stmt *ins = db_prep(v->db,
        "INSERT INTO ref_contig(tid,name,length,md5) VALUES(?,?,?,?);");
    int n = sam_hdr_nref(hdr);
    for (int i = 0; i < n; ++i) {
        const char *name = sam_hdr_tid2name(hdr, i);
        hts_pos_t   len  = sam_hdr_tid2len(hdr, i);
        kstring_t  m5    = { 0, 0, NULL };
        const char *md5  = NULL;
        if (sam_hdr_find_tag_pos((sam_hdr_t *)hdr, "SQ", i, "M5", &m5) == 0)
            md5 = m5.s;
        sqlite3_reset(ins);
        sqlite3_bind_int (ins, 1, i);
        sqlite3_bind_text(ins, 2, name, -1, SQLITE_STATIC);
        sqlite3_bind_int64(ins,3, (sqlite3_int64)len);
        if (md5) sqlite3_bind_text(ins, 4, md5, -1, SQLITE_TRANSIENT);
        else     sqlite3_bind_null(ins, 4);
        db_step_done(ins);
        free(m5.s);
    }
    sqlite3_finalize(ins);
}

sam_hdr_t *view_load_header(view_t *v) {
    if (v->backend_kind == VIEW_BACKEND_LMDB) {
#ifdef ISAM_ENABLE_LMDB
        char *txt = lmdb_backend_meta_get(v->lmdb, "alignment_header_text");
        if (!txt) die("view '%s' has no stored alignment header", v->path);
        sam_hdr_t *hdr = sam_hdr_init();
        if (!hdr) die("sam_hdr_init failed");
        if (sam_hdr_add_lines(hdr, txt, (int)strlen(txt)) < 0)
            die("sam_hdr_add_lines failed parsing stored header");
        free(txt);
        return hdr;
#else
        die("this build has no LMDB support");
#endif
    }

    sqlite3_stmt *s = db_prep(v->db,
        "SELECT value FROM meta WHERE key='alignment_header_text';");
    sam_hdr_t *hdr = NULL;
    if (db_step(s) == SQLITE_ROW) {
        const void *blob = sqlite3_column_blob(s, 0);
        int n = sqlite3_column_bytes(s, 0);
        hdr = sam_hdr_init();
        if (!hdr) die("sam_hdr_init failed");
        if (sam_hdr_add_lines(hdr, (const char *)blob, n) < 0)
            die("sam_hdr_add_lines failed parsing stored header");
    }
    sqlite3_finalize(s);
    if (!hdr)
        die("view '%s' has no stored alignment header", v->path);
    return hdr;
}

void view_link_store(view_t *v,
                     const char *store_path_override,
                     const char *store_index_path_override) {
    char *sp = store_path_override ? xstrdup(store_path_override)
                    : view_meta_get(v, "store_path");
    if (!sp) die("view '%s' has no store_path; pass --store", v->path);

    char *sip = store_index_path_override ? xstrdup(store_index_path_override)
                        : view_meta_get(v, "store_index_path");
    if (!sip) sip = isam_index_path_for(sp);

    v->store = store_open(sp, sip);
    free(sp); free(sip);
}

const char *view_tid2name(view_t *v, int tid) {
    if (!v->hdr || tid < 0) return NULL;
    return sam_hdr_tid2name(v->hdr, tid);
}

void view_meta_set(view_t *v, const char *key, const char *value) {
    if (v->backend_kind == VIEW_BACKEND_SQLITE) {
        meta_set(v->db, key, value);
        return;
    }
#ifdef ISAM_ENABLE_LMDB
    lmdb_view_put_meta(v->lmdb, key, value);
#else
    die("this build has no LMDB support");
#endif
}

void view_meta_set_i64(view_t *v, const char *key, int64_t value) {
    if (v->backend_kind == VIEW_BACKEND_SQLITE) {
        meta_set_i64(v->db, key, value);
        return;
    }
#ifdef ISAM_ENABLE_LMDB
    lmdb_view_put_meta_i64(v->lmdb, key, value);
#else
    die("this build has no LMDB support");
#endif
}

char *view_meta_get(view_t *v, const char *key) {
    if (v->backend_kind == VIEW_BACKEND_SQLITE)
        return meta_get(v->db, key);
#ifdef ISAM_ENABLE_LMDB
    return lmdb_backend_meta_get(v->lmdb, key);
#else
    die("this build has no LMDB support");
#endif
}

int view_meta_get_i64(view_t *v, const char *key, int64_t *out) {
    if (v->backend_kind == VIEW_BACKEND_SQLITE)
        return meta_get_i64(v->db, key, out);
#ifdef ISAM_ENABLE_LMDB
    return lmdb_backend_meta_get_i64(v->lmdb, key, out);
#else
    die("this build has no LMDB support");
#endif
}
