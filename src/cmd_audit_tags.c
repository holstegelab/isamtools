#include "common.h"
#include "aux_tags.h"

#include <htslib/hts.h>

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>

typedef struct audit_entry {
    uint16_t tag_code;
    int64_t count;
} audit_entry_t;

static void usage(void) {
    fputs(
        "Usage: isamtools audit-tags [options] <input.sam|bam|cram|->\n"
        "\n"
        "Inspect observed aux tags and classify them using built-in lists,\n"
        "header hints, and explicit tag-policy overrides.\n"
        "\n"
        "Options:\n"
        "  --share-tags TAGS          Tags reconstructed from canonical BAM/CRAM.\n"
        "  --no-share-tags TAGS       Comma-separated tags excluded from sharing.\n"
        "  --carry-tags TAGS          Comma-separated tags stored once per read.\n"
        "  --keep-alignment-tags TAGS Comma-separated tags kept per alignment.\n"
        "  -h, --help              Show this help.\n",
        stderr);
}

static audit_entry_t *audit_find_or_add(audit_entry_t **items,
                                        size_t *n,
                                        size_t *cap,
                                        const char tag[2]) {
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < *n; ++i) {
        if ((*items)[i].tag_code == code)
            return &(*items)[i];
    }
    if (*n == *cap) {
        size_t next_cap = *cap ? *cap * 2 : 16;
        *items = xrealloc(*items, next_cap * sizeof **items);
        *cap = next_cap;
    }
    (*items)[*n] = (audit_entry_t){ .tag_code = code, .count = 0 };
    return &(*items)[(*n)++];
}

static int audit_entry_cmp(const void *a, const void *b) {
    const audit_entry_t *x = a;
    const audit_entry_t *y = b;
    return x->tag_code < y->tag_code ? -1 : (x->tag_code > y->tag_code ? 1 : 0);
}

int cmd_audit_tags(int argc, char **argv) {
    const char *share_tags = NULL;
    const char *no_share_tags = NULL;
    const char *carry_tags = NULL;
    const char *keep_alignment_tags = NULL;

    static struct option longopts[] = {
        {"share-tags", required_argument, 0, 1001},
        {"no-share-tags", required_argument, 0, 1002},
        {"carry-tags", required_argument, 0, 1003},
        {"keep-alignment-tags", required_argument, 0, 1004},
        {"help",            no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1001: share_tags = optarg; break;
        case 1002: no_share_tags = optarg; break;
        case 1003: carry_tags = optarg; break;
        case 1004: keep_alignment_tags = optarg; break;
        case 'h': usage(); return 0;
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *input_path = argv[optind];

    samFile *input = sam_open(input_path, "r");
    if (!input) die("cannot open '%s'", input_path);
    sam_hdr_t *hdr = sam_hdr_read(input);
    if (!hdr) die("cannot read header from '%s'", input_path);

    isam_tag_classifier_t classifier;
    isam_tag_classifier_init(&classifier, hdr, NULL, share_tags, no_share_tags,
                             carry_tags, keep_alignment_tags);

    bam1_t *b = bam_init1();
    audit_entry_t *entries = NULL;
    size_t n_entries = 0, cap_entries = 0;
    int64_t n_records = 0;

    while (1) {
        int rc = sam_read1(input, hdr, b);
        if (rc == -1) break;
        if (rc < -1)
            die("error reading '%s' at record %lld", input_path, (long long)n_records);

        errno = 0;
        for (uint8_t *field = bam_aux_first(b); field; field = bam_aux_next(b, field)) {
            char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
            audit_entry_t *entry = audit_find_or_add(&entries, &n_entries, &cap_entries, tag);
            entry->count += 1;
        }
        if (errno && errno != ENOENT)
            die("corrupt aux data in '%s' at record %lld", input_path, (long long)n_records);
        ++n_records;
    }

    qsort(entries, n_entries, sizeof *entries, audit_entry_cmp);
    puts("tag\tscope\tsource\talignments");
    for (size_t i = 0; i < n_entries; ++i) {
        char tag[3];
        isam_tag_decode(entries[i].tag_code, tag);
        printf("%s\t%s\t%s\t%" PRId64 "\n",
               tag,
               isam_tag_scope_name(isam_tag_classifier_scope(&classifier, tag)),
               isam_tag_classifier_source(&classifier, tag),
               entries[i].count);
    }

    free(entries);
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(input);
    isam_tag_classifier_destroy(&classifier);
    return 0;
}
