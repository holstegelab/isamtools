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

static int has_path_suffix(const char *path, const char *suffix) {
    size_t path_len = strlen(path);
    size_t suffix_len = strlen(suffix);
    return path_len >= suffix_len &&
           strcmp(path + path_len - suffix_len, suffix) == 0;
}

char *isam_index_path_for(const char *store_path) {
    size_t stem_len = strlen(store_path);
    if (has_path_suffix(store_path, ".bam")) {
        stem_len -= 4;
    } else if (has_path_suffix(store_path, ".cram")) {
        stem_len -= 5;
    } else if (has_path_suffix(store_path, ".sam")) {
        stem_len -= 4;
    }

    char *r = xmalloc(stem_len + 6);
    memcpy(r, store_path, stem_len);
    memcpy(r + stem_len, ".isam", 6);
    return r;
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
