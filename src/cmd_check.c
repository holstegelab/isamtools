#include "common.h"
#include "aux_tags.h"
#include "store.h"
#include "view.h"
#include "rebuild.h"

#include <htslib/hts.h>
#include <getopt.h>
#include <errno.h>
#include <sys/stat.h>

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

static int input_aux_tags_are_present(const bam1_t *expected,
                                      const bam1_t *observed) {
    errno = 0;
    for (uint8_t *field = bam_aux_first(expected); field;
         field = bam_aux_next(expected, field)) {
        char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
        uint8_t *got = bam_aux_get(observed, tag);
        if (!got)
            return 0;
        int expected_len = isam_aux_field_payload_len(field);
        int got_len = isam_aux_field_payload_len(got);
        if (expected_len < 0 || got_len < 0)
            die("invalid aux tag '%s' while checking", tag);
        if (bam_aux_type(field) != bam_aux_type(got) ||
            expected_len != got_len ||
            memcmp(field + 1, got + 1, (size_t)expected_len) != 0)
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
        expected->core.qual != observed->core.qual ||
        expected->core.n_cigar != observed->core.n_cigar ||
        expected->core.mpos != observed->core.mpos ||
        expected->core.isize != observed->core.isize ||
        expected->core.l_qseq != observed->core.l_qseq)
        return 0;
    if (!same_ref_name(expected_hdr, expected->core.tid,
                       observed_hdr, observed->core.tid) ||
        !same_ref_name(expected_hdr, expected->core.mtid,
                       observed_hdr, observed->core.mtid))
        return 0;
    if (expected->core.n_cigar > 0 &&
        memcmp(bam_get_cigar(expected), bam_get_cigar(observed),
               (size_t)expected->core.n_cigar * sizeof(uint32_t)) != 0)
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

int isam_check_view_against_bam(const char *input_path,
                                const char *view_path,
                                const char *store_override,
                                const char *store_index_override) {
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
    bam1_t *input_b = bam_init1();
    bam1_t *canon_b = bam_init1();
    bam1_t *view_b = bam_init1();
    int64_t n_checked = 0;

    while (1) {
        int input_rc = sam_read1(input, input_hdr, input_b);
        lmdb_view_record_t rec;
        int view_has_record = lmdb_view_iter_next(it, &rec);
        if (input_rc == -1 && !view_has_record)
            break;
        if (input_rc == -1 && view_has_record)
            die("check: view '%s' has extra records after input ended", view_path);
        if (input_rc < -1)
            die("error reading input '%s' at record %lld",
                input_path, (long long)n_checked);
        if (!view_has_record)
            die("check: input '%s' has extra records after view ended", input_path);

        store_record_t srec;
        if (!store_lookup_read_id(v->store, rec.store_read_id, &srec))
            die("store_read_id '%lld' missing from store '%s'",
                (long long)rec.store_read_id, v->store->store_path);
        int frc = payload_cache_fetch(cache, v->store, &srec, canon_b);
        if (frc != 0)
            die("check: payload fetch failed for '%s' (rc=%d)", srec.qname, frc);

        isam_rebuild_bam(view_b, canon_b,
                         rec.flag, rec.tid, rec.pos0, rec.mapq,
                         rec.n_cigar, rec.cigar,
                         rec.mate_tid, rec.mate_pos0, rec.tlen,
                         rec.aux, rec.aux_len);
        lmdb_read_aux_apply(v->lmdb, rec.store_read_id, view_b, 0);
        lmdb_read_tags_apply(v->lmdb, rec.store_read_id, view_b, 0);
        isam_aux_tags_append_missing_from(canon_b, view_b,
                                          &shared_tags, &dropped_tags);
        lmdb_read_tags_apply(v->store->lmdb, rec.store_read_id, view_b, 0);
        isam_aux_tags_remove(view_b, &dropped_tags);

        if (!bam_record_matches_input(input_b, input_hdr, view_b, v->hdr))
            die("check: input record %lld (%s) differs from reconstructed view '%s'",
                (long long)(n_checked + 1), bam_get_qname(input_b), view_path);

        store_record_clear(&srec);
        lmdb_view_record_clear(&rec);
        ++n_checked;
    }

    bam_destroy1(input_b);
    bam_destroy1(canon_b);
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
    int allow_missing = 0;

    static struct option longopts[] = {
        {"view",          required_argument, 0, 1000},
        {"store",         required_argument, 0, 1001},
        {"store-index",   required_argument, 0, 1002},
        {"allow-missing", no_argument,       0, 1003},
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
        rc = isam_check_view_against_bam(input_path, view_path,
                                         store_path, store_index_path);
    } else {
        if (!store_path)
            die("--store is required when no alignment-view .isam is available");
        rc = isam_check_store_against_bam(input_path, store_path,
                                          store_index_path, allow_missing);
    }
    free(view_path_alloc);
    return rc;
}
