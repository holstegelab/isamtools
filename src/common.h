#ifndef ISAM_COMMON_H
#define ISAM_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

#define ISAM_VERSION       "0.1.0"
#define ISAM_SCHEMA_STORE  "1"
#define ISAM_SCHEMA_VIEW   "1"

void die(const char *fmt, ...) __attribute__((noreturn, format(printf, 1, 2)));
void warnx_(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

void *xmalloc(size_t n);
void *xcalloc(size_t n, size_t s);
void *xrealloc(void *p, size_t n);
char *xstrdup(const char *s);

/* qname hash: FNV-1a 64, returned as int64 (bit-cast). */
int64_t qname_hash(const char *qname);

int isam_has_path_suffix(const char *path, const char *suffix);
int isam_is_view_index_path(const char *path);
int isam_is_store_index_path(const char *path);

/* alignment view path helper: replace trailing .bam/.cram/.sam with .isam, else append .isam. caller frees. */
char *isam_index_path_for(const char *input_path);

/* canonical store index path helper: replace trailing .bam/.cram/.sam with .ssam, else append .ssam. caller frees. */
char *isam_store_index_path_for(const char *store_path);

/* iso8601 utc timestamp, malloc'd */
char *now_iso8601(void);

/* hash file header text (printable only). int64 fnv1a. */
int64_t header_text_hash(const char *text);

#endif
