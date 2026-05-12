#include "aux_tags.h"

#include <ctype.h>
#include <errno.h>

static const char *const builtin_read_tags[] = {
    "BC", "BZ", "CB", "CR", "CY", "FI", "FS", "FZ", "LB", "MI",
    "MM", "ML", "MN", "OX", "OQ", "PG", "PU", "QT", "QX", "RG",
    "RX"
};

static const char *const builtin_alignment_tags[] = {
    "AM", "AS", "BQ", "CC", "CG", "CM", "CP", "CQ", "CS", "E2",
    "H0", "H1", "H2", "HI", "IH", "MC", "MD", "MQ", "NH", "NM",
    "OA", "OC", "OP", "PQ", "Q2", "R2", "SA", "SM", "TC", "TS",
    "U2", "UQ", "tp", "cm", "s1", "s2", "ms", "nn", "dv", "de",
    "cg", "cs", "mc", "rm"
};

static const char *const builtin_carried_tags[] = {
    "HP", "PS"
};

static const char *const pacbio_read_tags[] = {
    "ac", "bc", "bl", "bq", "bt", "bx", "di", "dl", "dr", "ds",
    "ec", "ff", "fi", "fn", "fp", "gp", "ic", "ie", "il", "im",
    "ip", "is", "it", "iz", "lm", "ls", "lu", "ma", "nb", "nc",
    "np", "oc", "pw", "qe", "ql", "qs", "qt", "rc", "ri", "rn",
    "rp", "rq", "sa", "sm", "sn", "sx", "tm", "tu", "we", "ws",
    "zm"
};

static const char *const nanopore_read_tags[] = {
    "ch", "st", "du", "fn", "mx"
};

static int uint16_cmp(const void *a, const void *b) {
    uint16_t x = *(const uint16_t *)a;
    uint16_t y = *(const uint16_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static int tag_list_contains(const char *const *items, size_t n, const char tag[2]) {
    for (size_t i = 0; i < n; ++i) {
        if (items[i][0] == tag[0] && items[i][1] == tag[1])
            return 1;
    }
    return 0;
}

static void tag_set_add_static(isam_tag_set_t *set,
                               const char *const *items,
                               size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (isam_tag_set_add(set, items[i]) != 0)
            die("invalid built-in tag '%s'", items[i]);
    }
}

static int header_has_pacbio_platform(const sam_hdr_t *hdr) {
    const char *text = sam_hdr_str((sam_hdr_t *)hdr);
    if (!text) return 0;
    return strstr(text, "\tPL:PACBIO") != NULL || strstr(text, "\n@RG\tPL:PACBIO") != NULL;
}

static void header_add_manifest_tag(isam_tag_set_t *set,
                                    const char *text,
                                    const char *needle) {
    size_t needle_len = strlen(needle);
    const char *p = text;
    while ((p = strstr(p, needle)) != NULL) {
        const char *value = p + needle_len;
        if (value[0] && value[1] &&
            isam_tag_valid(value) &&
            (value[2] == 0 || value[2] == ';' || value[2] == '\t' || value[2] == '\n')) {
            char tag[3] = { value[0], value[1], 0 };
            if (isam_tag_set_add(set, tag) != 0)
                die("invalid tag '%s' in header manifest", tag);
        }
        p += needle_len;
    }
}

static void header_add_read_tags(isam_tag_classifier_t *classifier,
                                 const sam_hdr_t *hdr) {
    const char *text = sam_hdr_str((sam_hdr_t *)hdr);
    if (!text) return;

    classifier->pacbio_header = header_has_pacbio_platform(hdr);
    header_add_manifest_tag(&classifier->header_shared_tags, text, "Ipd:Frames=");
    header_add_manifest_tag(&classifier->header_shared_tags, text, "Ipd:CodecV1=");
    header_add_manifest_tag(&classifier->header_shared_tags, text, "PulseWidth:Frames=");
    header_add_manifest_tag(&classifier->header_shared_tags, text, "PulseWidth:CodecV1=");
}

static int aux_array_elem_size(char subtype) {
    switch (subtype) {
    case 'c': case 'C': return 1;
    case 's': case 'S': return 2;
    case 'i': case 'I': case 'f': return 4;
    case 'd': return 8;
    default: return -1;
    }
}

uint16_t isam_tag_code(const char tag[2]) {
    return (uint16_t)(((unsigned char)tag[0] << 8) | (unsigned char)tag[1]);
}

void isam_tag_decode(uint16_t code, char tag[3]) {
    tag[0] = (char)((code >> 8) & 0xff);
    tag[1] = (char)(code & 0xff);
    tag[2] = 0;
}

int isam_tag_valid(const char tag[2]) {
    return tag && isalnum((unsigned char)tag[0]) && isalnum((unsigned char)tag[1]);
}

void isam_tag_set_init(isam_tag_set_t *set) {
    memset(set, 0, sizeof *set);
}

void isam_tag_set_destroy(isam_tag_set_t *set) {
    free(set->items);
    memset(set, 0, sizeof *set);
}

int isam_tag_set_add(isam_tag_set_t *set, const char tag[2]) {
    if (!isam_tag_valid(tag)) return -1;
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < set->n; ++i) {
        if (set->items[i] == code) return 0;
    }
    if (set->n == set->cap) {
        size_t next_cap = set->cap ? set->cap * 2 : 16;
        set->items = xrealloc(set->items, next_cap * sizeof *set->items);
        set->cap = next_cap;
    }
    set->items[set->n++] = code;
    return 0;
}

int isam_tag_set_contains(const isam_tag_set_t *set, const char tag[2]) {
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < set->n; ++i) {
        if (set->items[i] == code) return 1;
    }
    return 0;
}

int isam_tag_set_parse_csv(isam_tag_set_t *set, const char *csv) {
    if (!csv || !csv[0]) return 0;
    const char *p = csv;
    while (*p) {
        while (*p == ',' || *p == ';' || isspace((unsigned char)*p)) ++p;
        if (!*p) break;
        const char *start = p;
        while (*p && *p != ',' && *p != ';' && !isspace((unsigned char)*p)) ++p;
        if ((size_t)(p - start) != 2) return -1;
        char tag[3] = { start[0], start[1], 0 };
        if (isam_tag_set_add(set, tag) != 0) return -1;
    }
    return 0;
}

char *isam_tag_set_join_csv(const isam_tag_set_t *set) {
    if (set->n == 0) return xstrdup("");
    uint16_t *copy = xmalloc(set->n * sizeof *copy);
    memcpy(copy, set->items, set->n * sizeof *copy);
    qsort(copy, set->n, sizeof *copy, uint16_cmp);

    size_t len = set->n * 2 + (set->n - 1) + 1;
    char *out = xmalloc(len);
    char *p = out;
    for (size_t i = 0; i < set->n; ++i) {
        char tag[3];
        isam_tag_decode(copy[i], tag);
        if (i) *p++ = ',';
        *p++ = tag[0];
        *p++ = tag[1];
    }
    *p = 0;
    free(copy);
    return out;
}

void isam_tag_classifier_init(isam_tag_classifier_t *classifier,
                              const sam_hdr_t *hdr,
                              const char *stored_shared_tags_csv,
                              const char *share_tags_csv,
                              const char *no_share_tags_csv,
                              const char *carry_tags_csv,
                              const char *keep_alignment_tags_csv) {
    memset(classifier, 0, sizeof *classifier);
    isam_tag_set_init(&classifier->configured_shared_tags);
    isam_tag_set_init(&classifier->no_share_tags);
    isam_tag_set_init(&classifier->header_shared_tags);
    isam_tag_set_init(&classifier->effective_shared_tags);
    isam_tag_set_init(&classifier->configured_carried_tags);
    isam_tag_set_init(&classifier->effective_carried_tags);
    isam_tag_set_init(&classifier->forced_alignment_tags);

    tag_set_add_static(&classifier->effective_shared_tags,
                       builtin_read_tags,
                       sizeof builtin_read_tags / sizeof builtin_read_tags[0]);
    tag_set_add_static(&classifier->effective_shared_tags,
                       nanopore_read_tags,
                       sizeof nanopore_read_tags / sizeof nanopore_read_tags[0]);
    tag_set_add_static(&classifier->effective_carried_tags,
                       builtin_carried_tags,
                       sizeof builtin_carried_tags / sizeof builtin_carried_tags[0]);
    if (stored_shared_tags_csv && isam_tag_set_parse_csv(&classifier->configured_shared_tags,
                                                         stored_shared_tags_csv) != 0)
        die("invalid stored shared tag list '%s'", stored_shared_tags_csv);
    if (share_tags_csv && isam_tag_set_parse_csv(&classifier->configured_shared_tags,
                                                 share_tags_csv) != 0)
        die("invalid shared tag list '%s'", share_tags_csv);
    if (no_share_tags_csv && isam_tag_set_parse_csv(&classifier->no_share_tags,
                                                    no_share_tags_csv) != 0)
        die("invalid no-share tag list '%s'", no_share_tags_csv);
    if (carry_tags_csv && isam_tag_set_parse_csv(&classifier->configured_carried_tags,
                                                 carry_tags_csv) != 0)
        die("invalid carried tag list '%s'", carry_tags_csv);
    if (keep_alignment_tags_csv && isam_tag_set_parse_csv(&classifier->forced_alignment_tags,
                                                          keep_alignment_tags_csv) != 0)
        die("invalid alignment-retained tag list '%s'", keep_alignment_tags_csv);

    if (hdr) header_add_read_tags(classifier, hdr);

    for (size_t i = 0; i < classifier->configured_shared_tags.n; ++i) {
        char tag[3];
        isam_tag_decode(classifier->configured_shared_tags.items[i], tag);
        isam_tag_set_add(&classifier->effective_shared_tags, tag);
    }
    for (size_t i = 0; i < classifier->header_shared_tags.n; ++i) {
        char tag[3];
        isam_tag_decode(classifier->header_shared_tags.items[i], tag);
        isam_tag_set_add(&classifier->effective_shared_tags, tag);
    }
    if (classifier->pacbio_header) {
        tag_set_add_static(&classifier->effective_shared_tags,
                           pacbio_read_tags,
                           sizeof pacbio_read_tags / sizeof pacbio_read_tags[0]);
    }
    for (size_t i = 0; i < classifier->configured_carried_tags.n; ++i) {
        char tag[3];
        isam_tag_decode(classifier->configured_carried_tags.items[i], tag);
        isam_tag_set_add(&classifier->effective_carried_tags, tag);
    }
    for (size_t i = 0; i < classifier->no_share_tags.n; ++i) {
        char tag[3];
        isam_tag_decode(classifier->no_share_tags.items[i], tag);
        for (size_t j = 0; j < classifier->effective_shared_tags.n; ++j) {
            if (classifier->effective_shared_tags.items[j] == classifier->no_share_tags.items[i]) {
                memmove(&classifier->effective_shared_tags.items[j],
                        &classifier->effective_shared_tags.items[j + 1],
                        (classifier->effective_shared_tags.n - j - 1) *
                            sizeof classifier->effective_shared_tags.items[0]);
                --classifier->effective_shared_tags.n;
                break;
            }
        }
        (void)tag;
    }
}

void isam_tag_classifier_destroy(isam_tag_classifier_t *classifier) {
    isam_tag_set_destroy(&classifier->configured_shared_tags);
    isam_tag_set_destroy(&classifier->no_share_tags);
    isam_tag_set_destroy(&classifier->header_shared_tags);
    isam_tag_set_destroy(&classifier->effective_shared_tags);
    isam_tag_set_destroy(&classifier->configured_carried_tags);
    isam_tag_set_destroy(&classifier->effective_carried_tags);
    isam_tag_set_destroy(&classifier->forced_alignment_tags);
    memset(classifier, 0, sizeof *classifier);
}

isam_tag_scope_t isam_tag_classifier_scope(const isam_tag_classifier_t *classifier,
                                           const char tag[2]) {
    if (isam_tag_set_contains(&classifier->forced_alignment_tags, tag))
        return ISAM_TAG_SCOPE_ALIGNMENT;
    if (isam_tag_set_contains(&classifier->effective_shared_tags, tag))
        return ISAM_TAG_SCOPE_SHARED;
    if (isam_tag_set_contains(&classifier->effective_carried_tags, tag))
        return ISAM_TAG_SCOPE_CARRIED;
    if (tag_list_contains(builtin_alignment_tags,
                          sizeof builtin_alignment_tags / sizeof builtin_alignment_tags[0],
                          tag))
        return ISAM_TAG_SCOPE_ALIGNMENT;
    return ISAM_TAG_SCOPE_UNKNOWN;
}

const char *isam_tag_scope_name(isam_tag_scope_t scope) {
    switch (scope) {
    case ISAM_TAG_SCOPE_SHARED: return "shared";
    case ISAM_TAG_SCOPE_CARRIED: return "carried";
    case ISAM_TAG_SCOPE_ALIGNMENT: return "alignment";
    default: return "unknown";
    }
}

const char *isam_tag_classifier_source(const isam_tag_classifier_t *classifier,
                                       const char tag[2]) {
    if (isam_tag_set_contains(&classifier->forced_alignment_tags, tag))
        return "configured-alignment";
    if (isam_tag_set_contains(&classifier->configured_shared_tags, tag))
        return "configured-shared";
    if (isam_tag_set_contains(&classifier->configured_carried_tags, tag))
        return "configured-carried";
    if (isam_tag_set_contains(&classifier->header_shared_tags, tag))
        return "header";
    if (tag_list_contains(builtin_read_tags,
                          sizeof builtin_read_tags / sizeof builtin_read_tags[0],
                          tag))
        return "builtin-shared";
    if (tag_list_contains(nanopore_read_tags,
                          sizeof nanopore_read_tags / sizeof nanopore_read_tags[0],
                          tag))
        return "builtin-shared";
    if (classifier->pacbio_header &&
        tag_list_contains(pacbio_read_tags,
                          sizeof pacbio_read_tags / sizeof pacbio_read_tags[0],
                          tag))
        return "pacbio";
    if (tag_list_contains(builtin_alignment_tags,
                          sizeof builtin_alignment_tags / sizeof builtin_alignment_tags[0],
                          tag))
        return "builtin-alignment";
    if (tag_list_contains(builtin_carried_tags,
                          sizeof builtin_carried_tags / sizeof builtin_carried_tags[0],
                          tag))
        return "builtin-carried";
    return "unknown";
}

int isam_aux_field_payload_len(const uint8_t *field) {
    char type = bam_aux_type(field);
    switch (type) {
    case 'A':
    case 'c':
    case 'C':
        return 1;
    case 's':
    case 'S':
        return 2;
    case 'i':
    case 'I':
    case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z':
    case 'H':
        return (int)strlen((const char *)(field + 1)) + 1;
    case 'B': {
        int elem_size = aux_array_elem_size((char)field[1]);
        if (elem_size < 0) {
            errno = EINVAL;
            return -1;
        }
        uint32_t len = bam_auxB_len(field);
        return 1 + 4 + (int)len * elem_size;
    }
    default:
        errno = EINVAL;
        return -1;
    }
}

void isam_aux_tags_append_missing_from(const bam1_t *src,
                                       bam1_t *dst,
                                       const isam_tag_set_t *include,
                                       const isam_tag_set_t *exclude) {
    if (!include || include->n == 0)
        return;
    errno = 0;
    for (uint8_t *field = bam_aux_first(src); field; field = bam_aux_next(src, field)) {
        char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
        if (!isam_tag_set_contains(include, tag))
            continue;
        if (exclude && isam_tag_set_contains(exclude, tag))
            continue;
        if (bam_aux_get(dst, tag))
            continue;
        int payload_len = isam_aux_field_payload_len(field);
        if (payload_len < 0)
            die("invalid canonical aux tag '%s'", tag);
        if (bam_aux_append(dst, tag, bam_aux_type(field), payload_len, field + 1) != 0)
            die("failed to append canonical aux tag '%s'", tag);
    }
    if (errno && errno != ENOENT)
        die("corrupt canonical aux data");
}

void isam_aux_tags_remove(bam1_t *record, const isam_tag_set_t *tags) {
    if (!tags || tags->n == 0)
        return;
    for (size_t i = 0; i < tags->n; ++i) {
        char tag[3];
        isam_tag_decode(tags->items[i], tag);
        uint8_t *field = bam_aux_get(record, tag);
        while (field) {
            if (bam_aux_del(record, field) != 0)
                die("failed to remove aux tag '%s'", tag);
            field = bam_aux_get(record, tag);
        }
    }
}
