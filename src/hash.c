#include "common.h"

int64_t qname_hash(const char *qname) {
    /* FNV-1a 64-bit -> bit-cast to int64 for SQLite INTEGER storage. */
    uint64_t h = 0xcbf29ce484222325ULL;
    for (const unsigned char *p = (const unsigned char *)qname; *p; ++p) {
        h ^= *p;
        h *= 0x100000001b3ULL;
    }
    union { uint64_t u; int64_t i; } cast = { .u = h };
    return cast.i;
}
