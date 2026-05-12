#ifndef ISAM_AUX_TAGS_H
#define ISAM_AUX_TAGS_H

#include "common.h"

typedef enum isam_tag_scope {
    ISAM_TAG_SCOPE_UNKNOWN = 0,
    ISAM_TAG_SCOPE_ALIGNMENT = 1,
    ISAM_TAG_SCOPE_CARRIED = 2,
    ISAM_TAG_SCOPE_SHARED = 3
} isam_tag_scope_t;

typedef struct isam_tag_set {
    uint16_t *items;
    size_t n;
    size_t cap;
} isam_tag_set_t;

typedef struct isam_tag_classifier {
    isam_tag_set_t configured_shared_tags;
    isam_tag_set_t no_share_tags;
    isam_tag_set_t header_shared_tags;
    isam_tag_set_t effective_shared_tags;
    isam_tag_set_t configured_carried_tags;
    isam_tag_set_t effective_carried_tags;
    isam_tag_set_t forced_alignment_tags;
    int pacbio_header;
} isam_tag_classifier_t;

uint16_t isam_tag_code(const char tag[2]);
void isam_tag_decode(uint16_t code, char tag[3]);
int isam_tag_valid(const char tag[2]);

void isam_tag_set_init(isam_tag_set_t *set);
void isam_tag_set_destroy(isam_tag_set_t *set);
int isam_tag_set_add(isam_tag_set_t *set, const char tag[2]);
int isam_tag_set_contains(const isam_tag_set_t *set, const char tag[2]);
int isam_tag_set_parse_csv(isam_tag_set_t *set, const char *csv);
char *isam_tag_set_join_csv(const isam_tag_set_t *set);

void isam_tag_classifier_init(isam_tag_classifier_t *classifier,
                              const sam_hdr_t *hdr,
                              const char *stored_shared_tags_csv,
                              const char *share_tags_csv,
                              const char *no_share_tags_csv,
                              const char *carry_tags_csv,
                              const char *keep_alignment_tags_csv);
void isam_tag_classifier_destroy(isam_tag_classifier_t *classifier);
isam_tag_scope_t isam_tag_classifier_scope(const isam_tag_classifier_t *classifier,
                                           const char tag[2]);
const char *isam_tag_scope_name(isam_tag_scope_t scope);
const char *isam_tag_classifier_source(const isam_tag_classifier_t *classifier,
                                       const char tag[2]);

int isam_aux_field_payload_len(const uint8_t *field);
void isam_aux_tags_append_missing_from(const bam1_t *src,
                                       bam1_t *dst,
                                       const isam_tag_set_t *include,
                                       const isam_tag_set_t *exclude);
void isam_aux_tags_remove(bam1_t *record, const isam_tag_set_t *tags);

#endif
