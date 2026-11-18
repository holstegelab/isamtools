#include "common.h"
#include "aux_tags.h"
#include "store.h"
#include "view.h"
#include "rebuild.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <errno.h>
#include <stdint.h>
#include <sys/stat.h>

#define DEFAULT_CHECK_BATCH_RECORDS ((size_t)65536)

static void usage(void) {
    fputs(
        "Usage: isamtools check [options] <input.sam|bam|cram>\n"
        "\n"
        "Verify a BAM/CRAM/SAM against an existing .isam view, or verify\n"
        "primary records against a canonical store.\n"
        "\n"
        "Options:\n"
        "  --view PATH         Alignment-view .isam to compare against.\n"
        "  --store PATH        Canonical store BAM/CRAM for store-only checks.\n"
        "  --store-index PATH  Store .ssam path.\n"
        "  --batch-records N   Group N view records by canonical CRAM container\n"
        "                      before fetching payloads (default 65536).\n"
        "  --allow-missing     Skip primary records absent from the store.\n"
        "  -h, --help          Show this help.\n",
        stderr);
}

static char *default_view_path_for(const char *input_path) {
    if (strcmp(input_path, "-") == 0)
        return NULL;
    return isam_index_path_for(input_path);
}

static int path_is_alignment_view(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0)
        return 0;
    if (S_ISDIR(st.st_mode))
        return 0;
    lmdb_backend_t *lmdb = lmdb_backend_open(path, 1, 0);
    char *type = lmdb_backend_meta_get(lmdb, "isam_type");
    int is_alignment = type && strcmp(type, "alignment") == 0;
    free(type);
    lmdb_backend_close(lmdb);
    return is_alignment;
}

static const char *record_ref_name(const sam_hdr_t *hdr, int tid) {
    if (tid < 0)
        return "*";
    const char *name = sam_hdr_tid2name((sam_hdr_t *)hdr, tid);
    return name ? name : "<missing>";
}

static int same_ref_name(const sam_hdr_t *a_hdr, int a_tid,
                         const sam_hdr_t *b_hdr, int b_tid) {
    return strcmp(record_ref_name(a_hdr, a_tid),
                  record_ref_name(b_hdr, b_tid)) == 0;
}

static int aux_type_is_integer(char type) {
    return type == 'c' || type == 'C' || type == 's' || type == 'S' ||
           type == 'i' || type == 'I';
}

static int aux_fields_are_equivalent(const uint8_t *expected,
                                     const uint8_t *observed) {
    char expected_type = bam_aux_type(expected);
    char observed_type = bam_aux_type(observed);
    if (aux_type_is_integer(expected_type) && aux_type_is_integer(observed_type))
        return bam_aux2i(expected) == bam_aux2i(observed);

    int expected_len = isam_aux_field_payload_len(expected);
    int observed_len = isam_aux_field_payload_len(observed);
    if (expected_len < 0 || observed_len < 0)
        return -1;
    return expected_type == observed_type && expected_len == observed_len &&
           memcmp(expected + 1, observed + 1, (size_t)expected_len) == 0;
}

static int cigar_op_is_matchlike(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

static uint64_t cigar_matchlike_block(const uint32_t *cigar, uint32_t n_cigar,
                                      uint32_t *offset, int *contains_m) {
    uint64_t length = 0;
    *contains_m = 0;
    while (*offset < n_cigar &&
           cigar_op_is_matchlike(bam_cigar_op(cigar[*offset]))) {
        if (bam_cigar_op(cigar[*offset]) == BAM_CMATCH)
            *contains_m = 1;
        length += bam_cigar_oplen(cigar[*offset]);
        ++*offset;
    }
    return length;
}

static int explicit_matchlike_blocks_are_equal(const uint32_t *a,
                                                uint32_t a_begin,
                                                uint32_t a_end,
                                                const uint32_t *b,
                                                uint32_t b_begin,
                                                uint32_t b_end) {
    uint32_t a_i = a_begin, b_i = b_begin;
    uint64_t a_left = 0, b_left = 0;
    uint32_t a_op = 0, b_op = 0;

    while (a_i < a_end || a_left > 0 || b_i < b_end || b_left > 0) {
        if (a_left == 0) {
            if (a_i == a_end)
                return 0;
            a_op = bam_cigar_op(a[a_i]);
            a_left = bam_cigar_oplen(a[a_i++]);
        }
        if (b_left == 0) {
            if (b_i == b_end)
                return 0;
            b_op = bam_cigar_op(b[b_i]);
            b_left = bam_cigar_oplen(b[b_i++]);
        }
        if (a_op != b_op)
            return 0;
        uint64_t consumed = a_left < b_left ? a_left : b_left;
        a_left -= consumed;
        b_left -= consumed;
    }
    return 1;
}

static int cigars_are_equivalent(const bam1_t *a, const bam1_t *b) {
    const uint32_t *a_cigar = bam_get_cigar(a);
    const uint32_t *b_cigar = bam_get_cigar(b);
    uint32_t a_i = 0, b_i = 0;

    while (a_i < a->core.n_cigar && b_i < b->core.n_cigar) {
        uint32_t a_op = bam_cigar_op(a_cigar[a_i]);
        uint32_t b_op = bam_cigar_op(b_cigar[b_i]);
        if (cigar_op_is_matchlike(a_op) || cigar_op_is_matchlike(b_op)) {
            if (!cigar_op_is_matchlike(a_op) || !cigar_op_is_matchlike(b_op))
                return 0;
            uint32_t a_begin = a_i, b_begin = b_i;
            int a_has_m, b_has_m;
            uint64_t a_len = cigar_matchlike_block(a_cigar, a->core.n_cigar,
                                                   &a_i, &a_has_m);
            uint64_t b_len = cigar_matchlike_block(b_cigar, b->core.n_cigar,
                                                   &b_i, &b_has_m);
            if (a_len != b_len)
                return 0;
            if (!a_has_m && !b_has_m &&
                !explicit_matchlike_blocks_are_equal(a_cigar, a_begin, a_i,
                                                     b_cigar, b_begin, b_i))
                return 0;
            continue;
        }

        if (a_op != b_op)
            return 0;
        uint64_t a_len = 0, b_len = 0;
        while (a_i < a->core.n_cigar &&
               bam_cigar_op(a_cigar[a_i]) == a_op)
            a_len += bam_cigar_oplen(a_cigar[a_i++]);
        while (b_i < b->core.n_cigar &&
               bam_cigar_op(b_cigar[b_i]) == b_op)
            b_len += bam_cigar_oplen(b_cigar[b_i++]);
        if (a_len != b_len)
            return 0;
    }
    return a_i == a->core.n_cigar && b_i == b->core.n_cigar;
}

static int mapping_qualities_are_equivalent(const bam1_t *a,
                                            const bam1_t *b) {
    if (a->core.qual == b->core.qual)
        return 1;

    /* HTSlib's CRAM decoder reports MAPQ 0 for unmapped records, even when
     * the input BAM used the SAM sentinel value 255 (MAPQ unavailable). */
    if ((a->core.flag & BAM_FUNMAP) && (b->core.flag & BAM_FUNMAP))
        return (a->core.qual == 255 && b->core.qual == 0) ||
               (a->core.qual == 0 && b->core.qual == 255);

    return 0;
}

static int input_aux_tags_are_present(const bam1_t *expected,
                                      const bam1_t *observed) {
    errno = 0;
    for (uint8_t *field = bam_aux_first(expected); field;
         field = bam_aux_next(expected, field)) {
        char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
        uint8_t *got = bam_aux_get(observed, tag);
        if (!got)
            return 0;
        int equivalent = aux_fields_are_equivalent(field, got);
        if (equivalent < 0)
            die("invalid aux tag '%s' while checking", tag);
        if (!equivalent)
            return 0;
    }
    if (errno && errno != ENOENT)
        die("corrupt aux data while checking");
    return 1;
}

static int bam_record_matches_input(const bam1_t *expected,
                                    const sam_hdr_t *expected_hdr,
                                    const bam1_t *observed,
                                    const sam_hdr_t *observed_hdr) {
    if (strcmp(bam_get_qname(expected), bam_get_qname(observed)) != 0)
        return 0;
    if ((expected->core.flag & 0xffff) != (observed->core.flag & 0xffff) ||
        expected->core.pos != observed->core.pos ||
        !mapping_qualities_are_equivalent(expected, observed) ||
        expected->core.mpos != observed->core.mpos ||
        expected->core.isize != observed->core.isize ||
        expected->core.l_qseq != observed->core.l_qseq)
        return 0;
    if (!same_ref_name(expected_hdr, expected->core.tid,
                       observed_hdr, observed->core.tid) ||
        !same_ref_name(expected_hdr, expected->core.mtid,
                       observed_hdr, observed->core.mtid))
        return 0;
    if (!cigars_are_equivalent(expected, observed))
        return 0;
    size_t seq_bytes = (size_t)((expected->core.l_qseq + 1) / 2);
    if (seq_bytes > 0 &&
        memcmp(bam_get_seq(expected), bam_get_seq(observed), seq_bytes) != 0)
        return 0;
    if (expected->core.l_qseq > 0 &&
        memcmp(bam_get_qual(expected), bam_get_qual(observed),
               (size_t)expected->core.l_qseq) != 0)
        return 0;
    return input_aux_tags_are_present(expected, observed);
}

static void load_view_tag_sets(view_t *v,
                               isam_tag_set_t *shared_tags,
                               isam_tag_set_t *dropped_tags) {
    isam_tag_set_init(shared_tags);
    isam_tag_set_init(dropped_tags);
    char *shared_csv = view_meta_get(v, "shared_tags_snapshot");
    char *dropped_csv = view_meta_get(v, "dropped_tags");
    if (isam_tag_set_parse_csv(shared_tags, shared_csv) != 0)
        die("invalid shared tag metadata in '%s'", v->path);
    if (isam_tag_set_parse_csv(dropped_tags, dropped_csv) != 0)
        die("invalid dropped tag metadata in '%s'", v->path);
    free(shared_csv);
    free(dropped_csv);
}

typedef struct check_batch_record {
    bam1_t *input;
    bam1_t *canonical;
    lmdb_view_record_t view;
    store_record_t store;
} check_batch_record_t;

typedef struct check_fetch_order {
    size_t batch_index;
    int has_cram_loc;
    int64_t cram_coffset;
    int32_t cram_within;
    int64_t read_id;
} check_fetch_order_t;

static int compare_check_fetch_order(const void *va, const void *vb) {
    const check_fetch_order_t *a = va;
    const check_fetch_order_t *b = vb;
    if (a->has_cram_loc != b->has_cram_loc)
        return a->has_cram_loc ? -1 : 1;
    if (a->has_cram_loc) {
        if (a->cram_coffset != b->cram_coffset)
            return a->cram_coffset < b->cram_coffset ? -1 : 1;
        if (a->cram_within != b->cram_within)
            return a->cram_within < b->cram_within ? -1 : 1;
        if (a->read_id != b->read_id)
            return a->read_id < b->read_id ? -1 : 1;
    }
    if (a->batch_index == b->batch_index)
        return 0;
    return a->batch_index < b->batch_index ? -1 : 1;
}

static size_t parse_positive_size(const char *text, const char *what) {
    char *end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno || !end || *end || value == 0 || value > SIZE_MAX)
        die("invalid %s '%s' (expected a positive integer)", what, text);
    return (size_t)value;
}

static int isam_check_view_against_bam_batched(
                                const char *input_path,
                                const char *view_path,
                                const char *store_override,
                                const char *store_index_override,
                                size_t batch_records) {
    if (batch_records == 0)
        die("check batch size must be positive");
    if (batch_records > SIZE_MAX / sizeof(check_batch_record_t) ||
        batch_records > SIZE_MAX / sizeof(check_fetch_order_t))
        die("check batch size is too large");
    samFile *input = sam_open(input_path, "r");
    if (!input)
        die("cannot open input '%s'", input_path);
    sam_hdr_t *input_hdr = sam_hdr_read(input);
    if (!input_hdr)
        die("cannot read header from '%s'", input_path);

    view_t *v = view_open_read(view_path);
    view_link_store(v, store_override, store_index_override);
    isam_tag_set_t shared_tags, dropped_tags;
    load_view_tag_sets(v, &shared_tags, &dropped_tags);

    lmdb_view_iter_t *it = lmdb_view_iter_open(v->lmdb, LMDB_VIEW_ITER_IMPORT,
                                               -1, 0, 0);
    payload_cache_t *cache = payload_cache_create(64 * 1024 * 1024);
    bam1_t *view_b = bam_init1();
    check_batch_record_t *batch = xcalloc(batch_records, sizeof *batch);
    check_fetch_order_t *fetch_order =
        xmalloc(batch_records * sizeof *fetch_order);
    int64_t n_checked = 0;

    while (1) {
        size_t n_batch = 0;
        int reached_end = 0;
        while (n_batch < batch_records) {
            check_batch_record_t *item = &batch[n_batch];
            if (!item->input) {
                item->input = bam_init1();
                item->canonical = bam_init1();
                if (!item->input || !item->canonical)
                    die("bam_init1 failed");
            }
            int input_rc = sam_read1(input, input_hdr, item->input);
            int view_has_record = lmdb_view_iter_next(it, &item->view);
            if (input_rc == -1 && !view_has_record) {
                reached_end = 1;
                break;
            }
            if (input_rc == -1 && view_has_record)
                die("check: view '%s' has extra records after input ended", view_path);
            if (input_rc < -1)
                die("error reading input '%s' at record %lld",
                    input_path, (long long)(n_checked + (int64_t)n_batch));
            if (!view_has_record)
                die("check: input '%s' has extra records after view ended", input_path);

            if (!store_lookup_read_id(v->store, item->view.store_read_id,
                                      &item->store))
                die("store_read_id '%lld' missing from store '%s'",
                    (long long)item->view.store_read_id,
                    v->store->store_path);
            fetch_order[n_batch] = (check_fetch_order_t) {
                .batch_index = n_batch,
                .has_cram_loc = item->store.has_cram_loc,
                .cram_coffset = item->store.cram_coffset,
                .cram_within = item->store.cram_within,
                .read_id = item->store.read_id,
            };
            ++n_batch;
        }

        if (n_batch == 0)
            break;

        qsort(fetch_order, n_batch, sizeof *fetch_order,
              compare_check_fetch_order);
        for (size_t i = 0; i < n_batch; ++i) {
            check_batch_record_t *item = &batch[fetch_order[i].batch_index];
            int frc = payload_cache_fetch(cache, v->store, &item->store,
                                          item->canonical);
            if (frc != 0)
                die("check: payload fetch failed for '%s' (rc=%d)",
                    item->store.qname, frc);
        }

        for (size_t i = 0; i < n_batch; ++i) {
            check_batch_record_t *item = &batch[i];
            lmdb_view_record_t *rec = &item->view;
            isam_rebuild_bam(view_b, item->canonical,
                             rec->flag, rec->tid, rec->pos0, rec->mapq,
                             rec->n_cigar, rec->cigar,
                             rec->mate_tid, rec->mate_pos0, rec->tlen,
                             rec->aux, rec->aux_len);
            lmdb_read_aux_apply(v->lmdb, rec->store_read_id, view_b, 0);
            lmdb_read_tags_apply(v->lmdb, rec->store_read_id, view_b, 0);
            isam_aux_tags_append_missing_from(item->canonical, view_b,
                                              &shared_tags, &dropped_tags);
            lmdb_read_tags_apply(v->store->lmdb, rec->store_read_id,
                                 view_b, 0);
            isam_aux_tags_remove(view_b, &dropped_tags);

            if (!bam_record_matches_input(item->input, input_hdr,
                                          view_b, v->hdr))
                die("check: input record %lld (%s) differs from reconstructed view '%s'",
                    (long long)(n_checked + (int64_t)i + 1),
                    bam_get_qname(item->input), view_path);

            store_record_clear(&item->store);
            lmdb_view_record_clear(rec);
        }
        n_checked += (int64_t)n_batch;
        if (reached_end)
            break;
    }

    for (size_t i = 0; i < batch_records; ++i) {
        if (batch[i].input)
            bam_destroy1(batch[i].input);
        if (batch[i].canonical)
            bam_destroy1(batch[i].canonical);
    }
    free(fetch_order);
    free(batch);
    bam_destroy1(view_b);
    payload_cache_destroy(cache);
    lmdb_view_iter_close(it);
    isam_tag_set_destroy(&shared_tags);
    isam_tag_set_destroy(&dropped_tags);
    view_close(v);
    sam_hdr_destroy(input_hdr);
    sam_close(input);

    fprintf(stderr,
            "isamtools check: %lld records match reconstructed view '%s'\n",
            (long long)n_checked, view_path);
    return 0;
}

int isam_check_view_against_bam(const char *input_path,
                                const char *view_path,
                                const char *store_override,
                                const char *store_index_override) {
    return isam_check_view_against_bam_batched(
        input_path, view_path, store_override, store_index_override,
        DEFAULT_CHECK_BATCH_RECORDS);
}

int isam_check_store_against_bam(const char *input_path,
                                 const char *store_path,
                                 const char *store_index_path,
                                 int allow_missing) {
    char *store_idx_alloc = NULL;
    if (!store_index_path) {
        store_idx_alloc = isam_store_index_path_for(store_path);
        store_index_path = store_idx_alloc;
    }
    store_t *store = store_open(store_path, store_index_path);
    if (!store_has_name_index(store))
        die("store index '%s' has no QNAME index", store_index_path);

    samFile *input = sam_open(input_path, "r");
    if (!input)
        die("cannot open input '%s'", input_path);
    sam_hdr_t *input_hdr = sam_hdr_read(input);
    if (!input_hdr)
        die("cannot read header from '%s'", input_path);

    payload_cache_t *cache = payload_cache_create(64 * 1024 * 1024);
    bam1_t *input_b = bam_init1();
    bam1_t *store_b = bam_init1();
    int64_t n_records = 0, n_checked = 0, n_skipped = 0, n_missing = 0;

    while (1) {
        int rc = sam_read1(input, input_hdr, input_b);
        if (rc == -1) break;
        if (rc < -1)
            die("error reading input '%s' at record %lld",
                input_path, (long long)n_records);
        ++n_records;

        if (input_b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            ++n_skipped;
            continue;
        }

        const char *qname = bam_get_qname(input_b);
        store_record_t srec;
        if (!store_lookup_qname(store, qname, &srec)) {
            if (!allow_missing)
                die("check: QNAME '%s' not found in canonical store index '%s'",
                    qname, store_path);
            ++n_missing;
            continue;
        }

        int frc = payload_cache_fetch(cache, store, &srec, store_b);
        if (frc != 0)
            die("check: payload fetch failed for '%s' (rc=%d)", qname, frc);
        if (!bam_record_matches_input(input_b, input_hdr, store_b, store->hdr))
            die("check: primary record for QNAME '%s' differs from store '%s'",
                qname, store_path);

        store_record_clear(&srec);
        ++n_checked;
    }

    bam_destroy1(input_b);
    bam_destroy1(store_b);
    payload_cache_destroy(cache);
    sam_hdr_destroy(input_hdr);
    sam_close(input);
    store_close(store);
    free(store_idx_alloc);

    fprintf(stderr,
            "isamtools check: %lld primary records match store, %lld secondary/supplementary skipped, %lld missing skipped\n",
            (long long)n_checked,
            (long long)n_skipped,
            (long long)n_missing);
    return 0;
}

int cmd_check(int argc, char **argv) {
    const char *view_path = NULL;
    const char *store_path = NULL;
    const char *store_index_path = NULL;
    size_t batch_records = DEFAULT_CHECK_BATCH_RECORDS;
    int batch_records_set = 0;
    int allow_missing = 0;

    static struct option longopts[] = {
        {"view",          required_argument, 0, 1000},
        {"store",         required_argument, 0, 1001},
        {"store-index",   required_argument, 0, 1002},
        {"allow-missing", no_argument,       0, 1003},
        {"batch-records", required_argument, 0, 1004},
        {"help",          no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: view_path = optarg; break;
        case 1001: store_path = optarg; break;
        case 1002: store_index_path = optarg; break;
        case 1003: allow_missing = 1; break;
        case 1004:
            batch_records = parse_positive_size(optarg, "--batch-records");
            batch_records_set = 1;
            break;
        case 'h': usage(); return 0;
        default: usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *input_path = argv[optind];

    char *view_path_alloc = NULL;
    if (!view_path) {
        view_path_alloc = default_view_path_for(input_path);
        if (view_path_alloc && path_is_alignment_view(view_path_alloc))
            view_path = view_path_alloc;
    }

    int rc;
    if (view_path) {
        rc = isam_check_view_against_bam_batched(
            input_path, view_path, store_path, store_index_path,
            batch_records);
    } else {
        if (batch_records_set)
            die("--batch-records requires an alignment-view check");
        if (!store_path)
            die("--store is required when no alignment-view .isam is available");
        rc = isam_check_store_against_bam(input_path, store_path,
                                          store_index_path, allow_missing);
    }
    free(view_path_alloc);
    return rc;
}
