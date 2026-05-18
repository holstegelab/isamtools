#include "view.h"

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
    v->lmdb = lmdb_backend_open(path, /*readonly*/1, 0);

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
    view_t *v = xcalloc(1, sizeof *v);
    v->path = xstrdup(path);
    struct stat st;
    if (stat(path, &st) == 0)
        die("LMDB output path '%s' already exists; choose a new path", path);
    if (errno != ENOENT)
        die("stat '%s': %s", path, strerror(errno));
    v->lmdb = lmdb_backend_open(path, /*readonly*/0, 0);
    return v;
}

void view_close(view_t *v) {
    if (!v) return;
    if (v->store) store_close(v->store);
    if (v->hdr)   sam_hdr_destroy(v->hdr);
    if (v->lmdb)  lmdb_backend_close(v->lmdb);
    free(v->path);
    free(v);
}

void view_set_header(view_t *v, const sam_hdr_t *hdr) {
    /* persist header text */
    const char *txt = sam_hdr_str((sam_hdr_t *)hdr);
    if (!txt) die("sam_hdr_str returned NULL");
    lmdb_view_put_meta(v->lmdb, "alignment_header_text", txt);
    lmdb_view_put_meta_i64(v->lmdb, "alignment_header_hash",
                           header_text_hash(txt));
}

sam_hdr_t *view_load_header(view_t *v) {
    char *txt = lmdb_backend_meta_get(v->lmdb, "alignment_header_text");
    if (!txt) die("view '%s' has no stored alignment header", v->path);
    sam_hdr_t *hdr = sam_hdr_init();
    if (!hdr) die("sam_hdr_init failed");
    if (sam_hdr_add_lines(hdr, txt, (int)strlen(txt)) < 0)
        die("sam_hdr_add_lines failed parsing stored header");
    free(txt);
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
    if (!sip) sip = isam_store_index_path_for(sp);

    v->store = store_open(sp, sip);
    free(sp); free(sip);
}

const char *view_tid2name(view_t *v, int tid) {
    if (!v->hdr || tid < 0) return NULL;
    return sam_hdr_tid2name(v->hdr, tid);
}

void view_meta_set(view_t *v, const char *key, const char *value) {
    lmdb_view_put_meta(v->lmdb, key, value);
}

void view_meta_set_i64(view_t *v, const char *key, int64_t value) {
    lmdb_view_put_meta_i64(v->lmdb, key, value);
}

char *view_meta_get(view_t *v, const char *key) {
    return lmdb_backend_meta_get(v->lmdb, key);
}

int view_meta_get_i64(view_t *v, const char *key, int64_t *out) {
    return lmdb_backend_meta_get_i64(v->lmdb, key, out);
}
