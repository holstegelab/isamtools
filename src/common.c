#include <stdarg.h>
#include <time.h>

#include "common.h"

void die(const char *fmt, ...) {
    fputs("isamtools: ", stderr);
    va_list ap; va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

void warnx_(const char *fmt, ...) {
    fputs("isamtools: ", stderr);
    va_list ap; va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}
void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n ? n : 1, s ? s : 1);
    if (!p) die("out of memory");
    return p;
}
void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) die("out of memory");
    return q;
}
char *xstrdup(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *r = xmalloc(n);
    memcpy(r, s, n);
    return r;
}

int isam_has_path_suffix(const char *path, const char *suffix) {
    size_t path_len = strlen(path);
    size_t suffix_len = strlen(suffix);
    return path_len >= suffix_len &&
           strcmp(path + path_len - suffix_len, suffix) == 0;
}

int isam_is_view_index_path(const char *path) {
    return isam_has_path_suffix(path, ".isam");
}

int isam_is_store_index_path(const char *path) {
    return isam_has_path_suffix(path, ".ssam");
}

static char *index_path_for(const char *input_path, const char *extension) {
    size_t stem_len = strlen(input_path);
    if (isam_has_path_suffix(input_path, ".bam")) {
        stem_len -= 4;
    } else if (isam_has_path_suffix(input_path, ".cram")) {
        stem_len -= 5;
    } else if (isam_has_path_suffix(input_path, ".sam")) {
        stem_len -= 4;
    }

    size_t ext_len = strlen(extension);
    char *r = xmalloc(stem_len + ext_len + 1);
    memcpy(r, input_path, stem_len);
    memcpy(r + stem_len, extension, ext_len + 1);
    return r;
}

char *isam_index_path_for(const char *input_path) {
    return index_path_for(input_path, ".isam");
}

char *isam_store_index_path_for(const char *store_path) {
    return index_path_for(store_path, ".ssam");
}

char *now_iso8601(void) {
    char buf[40];
    time_t t = time(NULL);
    struct tm tm;
    gmtime_r(&t, &tm);
    strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
    return xstrdup(buf);
}

int64_t header_text_hash(const char *text) {
    if (!text) return 0;
    /* FNV-1a 64-bit */
    uint64_t h = 0xcbf29ce484222325ULL;
    for (const unsigned char *p = (const unsigned char *)text; *p; ++p) {
        h ^= *p;
        h *= 0x100000001b3ULL;
    }
    union { uint64_t u; int64_t i; } cast = { .u = h };
    return cast.i;
}
