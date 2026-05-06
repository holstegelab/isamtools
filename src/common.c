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

char *isam_index_path_for(const char *store_path) {
    size_t n = strlen(store_path);
    char *r = xmalloc(n + 6);
    memcpy(r, store_path, n);
    memcpy(r + n, ".isam", 6);
    return r;
}

int isam_file_is_sqlite(const char *path) {
    unsigned char magic[16];
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    size_t n = fread(magic, 1, sizeof magic, f);
    fclose(f);
    static const unsigned char sqlite_magic[16] = {
        'S','Q','L','i','t','e',' ','f','o','r','m','a','t',' ','3','\0'
    };
    return n == sizeof magic && memcmp(magic, sqlite_magic, sizeof magic) == 0;
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
