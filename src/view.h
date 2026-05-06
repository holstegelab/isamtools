#ifndef ISAM_VIEW_H
#define ISAM_VIEW_H

#include "common.h"
#include "store.h"

typedef struct view {
    char       *path;
    lmdb_backend_t *lmdb;
    sam_hdr_t  *hdr;            /* alignment-view header (from meta) */
    int         n_targets;
    /* tid->name not stored; we use sam_hdr_tid2name via hdr */
    store_t    *store;          /* may be NULL until linked */
} view_t;

/* Open an existing view file for reading (db read-only). hdr is loaded. */
view_t *view_open_read(const char *path);

/* Open or create a view file for writing. Applies VIEW_DDL on create. */
view_t *view_open_write(const char *path);

void    view_close(view_t *v);

/* Save SAM header (text form) into meta as 'alignment_header_text', and
 * populate ref_contig from header SQ lines. */
void    view_set_header(view_t *v, const sam_hdr_t *hdr);

void    view_meta_set(view_t *v, const char *key, const char *value);
void    view_meta_set_i64(view_t *v, const char *key, int64_t value);
char   *view_meta_get(view_t *v, const char *key);
int     view_meta_get_i64(view_t *v, const char *key, int64_t *out);

/* Read the SAM header from meta. Caller does not own; freed with view. */
sam_hdr_t *view_load_header(view_t *v);

/* Resolve / open the canonical store referenced by this view's metadata,
 * with optional overrides (NULL = use stored values). Sets v->store. */
void    view_link_store(view_t *v,
                        const char *store_path_override,
                        const char *store_index_path_override);

/* Translate tid in this view's header into a contig name. */
const char *view_tid2name(view_t *v, int tid);

#endif
