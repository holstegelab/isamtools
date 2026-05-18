#include "common.h"
#include "aux_tags.h"
#include "store.h"
#include "view.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <getopt.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

int isam_check_view_against_bam(const char *input_path,
                                const char *view_path,
                                const char *store_override,
                                const char *store_index_override);

static void usage(void) {
    fputs(
           "Usage: isamtools import [options] <input.sam|bam|cram|-> \\\n"
        "         --store <store.bam|cram> [--aln <out.isam>]\n"
        "\n"
        "Options:\n"
        "  --store PATH         Canonical store (BAM/CRAM).\n"
        "  --store-index PATH   Store .isam (default <store-root>.isam).\n"
        "  --aln PATH           Output alignment-view .isam (default <input-root>.isam).\n"
        "  --name NAME          Logical alias for the alignment view.\n"
        "  --carry-tags TAGS    Comma-separated source tags stored once per read.\n"
        "  --keep-alignment-tags TAGS  Comma-separated source tags kept inline per alignment.\n"
        "  --drop-tags TAGS     Suppress tags from source and canonical/store reconstruction.\n"
        "  --drop-source-tags TAGS  Ignore matching tags from the imported source only.\n"
        "  --check            After import, verify input records match the reconstructed view.\n"
        "  --strict-tag-policy  Fail on stored-once tag conflicts.\n"
        "  --verbose-tag-policy Print a detailed tag policy summary.\n"
        "  --allow-missing      Skip records whose QNAME isn't in the store.\n"
        "  --write-index        Run `sort --all` after import.\n"
        "  -h, --help           Show this help.\n",
        stderr);
}

KHASH_SET_INIT_INT64(u64set)

typedef struct tag_count {
    uint16_t tag_code;
    int64_t count;
} tag_count_t;

typedef struct read_aux_dict_trainer {
    kstring_t samples;
    size_t sample_sizes[1000];
    unsigned sample_count;
    int trained;
} read_aux_dict_trainer_t;

static void tag_count_inc(tag_count_t **items, size_t *n, size_t *cap,
                          const char tag[2]) {
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < *n; ++i) {
        if ((*items)[i].tag_code == code) {
            ++(*items)[i].count;
            return;
        }
    }
    if (*n == *cap) {
        size_t next_cap = *cap ? *cap * 2 : 8;
        *items = xrealloc(*items, next_cap * sizeof **items);
        *cap = next_cap;
    }
    (*items)[*n] = (tag_count_t){ .tag_code = code, .count = 1 };
    ++*n;
}

static uint64_t read_tag_seen_key(int64_t read_id, const char tag[2]) {
    return ((uint64_t)read_id << 16) | (uint64_t)isam_tag_code(tag);
}

static int read_tag_seen_add(khash_t(u64set) *set,
                             int64_t read_id,
                             const char tag[2]) {
    int absent = 0;
    khiter_t iter = kh_put(u64set, set, read_tag_seen_key(read_id, tag), &absent);
    if (absent < 0)
        die("out of memory");
    (void)iter;
    return absent;
}

static int read_tag_seen_contains(khash_t(u64set) *set,
                                  int64_t read_id,
                                  const char tag[2]) {
    return kh_get(u64set, set, read_tag_seen_key(read_id, tag)) != kh_end(set);
}

static int aux_field_value_matches(const uint8_t *field,
                                   char type,
                                   const void *value,
                                   size_t value_len) {
    if (!field)
        return 0;
    int payload_len = isam_aux_field_payload_len(field);
    if (payload_len < 0)
        die("invalid canonical aux tag while comparing import tags");
    return bam_aux_type(field) == type &&
           (size_t)payload_len == value_len &&
           (value_len == 0 || memcmp(field + 1, value, value_len) == 0);
}

static const uint8_t *aux_blob_find_field(const uint8_t *aux,
                                          size_t aux_len,
                                          const char tag[2]) {
    const uint8_t *p = aux;
    const uint8_t *end = aux + aux_len;
    while (p < end) {
        if (p + 3 > end)
            die("truncated read aux blob");
        const uint8_t *field = p + 2;
        int payload_len = isam_aux_field_payload_len(field);
        if (payload_len < 0 || p + 3 + payload_len > end)
            die("corrupt read aux blob");
        if (p[0] == (uint8_t)tag[0] && p[1] == (uint8_t)tag[1])
            return field;
        p += 3 + payload_len;
    }
    return NULL;
}

static void aux_blob_append_field(uint8_t **aux,
                                  size_t *aux_len,
                                  size_t *aux_cap,
                                  const uint8_t *field) {
    int payload_len = isam_aux_field_payload_len(field);
    if (payload_len < 0)
        die("invalid aux tag while bundling read aux");
    size_t field_len = 3 + (size_t)payload_len;
    if (*aux_len + field_len > *aux_cap) {
        size_t next = *aux_cap ? *aux_cap * 2 : 64;
        while (next < *aux_len + field_len)
            next *= 2;
        *aux = xrealloc(*aux, next);
        *aux_cap = next;
    }
    memcpy(*aux + *aux_len, field - 2, field_len);
    *aux_len += field_len;
}

static void maybe_train_read_aux_dict(lmdb_backend_t *lmdb,
                                      read_aux_dict_trainer_t *trainer,
                                      const uint8_t *aux,
                                      size_t aux_len) {
    if (trainer->trained || trainer->sample_count >= 1000 || aux_len == 0)
        return;
    if (aux_len > INT_MAX)
        die("read aux blob is too large for dictionary training");
    kputsn((const char *)aux, (int)aux_len, &trainer->samples);
    trainer->sample_sizes[trainer->sample_count++] = aux_len;
    if (trainer->sample_count == 1000) {
        trainer->trained = lmdb_read_aux_set_zstd_dict(lmdb,
                                                       trainer->samples.s,
                                                       trainer->sample_sizes,
                                                       trainer->sample_count,
                                                       64 * 1024);
        if (trainer->trained)
            fprintf(stderr,
                    "isamtools import: trained read_aux zstd dictionary from 1000 records (%zu bytes)\n",
                    lmdb->read_aux_zstd_dict_len);
    }
}

static char *lock_path_for(const char *path) {
    size_t n = strlen(path);
    char *lock_path = xmalloc(n + 6);
    memcpy(lock_path, path, n);
    memcpy(lock_path + n, "-lock", 6);
    return lock_path;
}

static char *default_aln_path_for(const char *input_path) {
    if (strcmp(input_path, "-") == 0)
        return NULL;
    return isam_index_path_for(input_path);
}

static int import_output_is_incomplete(const char *path) {
    lmdb_backend_t *lmdb = lmdb_backend_open(path, /*readonly*/1, 0);
    char *type = lmdb_backend_meta_get(lmdb, "isam_type");
    if (!type || strcmp(type, "alignment") != 0) {
        free(type);
        lmdb_backend_close(lmdb);
        return 0;
    }

    char *state = lmdb_backend_meta_get(lmdb, "import_state");
    int incomplete = 0;
    if (state && strcmp(state, "in_progress") == 0) {
        incomplete = 1;
    } else {
        char *max_ref_span = lmdb_backend_meta_get(lmdb, "view_max_ref_span");
        char *max_query_len = lmdb_backend_meta_get(lmdb, "view_max_query_len");
        char *tag_policy_version = lmdb_backend_meta_get(lmdb, "tag_policy_version");
        incomplete = !max_ref_span || !max_query_len || !tag_policy_version;
        free(max_ref_span);
        free(max_query_len);
        free(tag_policy_version);
    }

    free(state);
    free(type);
    lmdb_backend_close(lmdb);
    return incomplete;
}

static void prepare_import_output_path(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) {
        if (errno == ENOENT)
            return;
        die("stat '%s': %s", path, strerror(errno));
    }

    if (S_ISDIR(st.st_mode))
        die("LMDB output path '%s' is a directory", path);

    if (!import_output_is_incomplete(path))
        die("LMDB output path '%s' already exists; choose a new path", path);

    warnx_("import: removing incomplete output '%s' from previous failed import", path);
    if (unlink(path) != 0 && errno != ENOENT)
        die("unlink '%s': %s", path, strerror(errno));

    char *lock_path = lock_path_for(path);
    if (unlink(lock_path) != 0 && errno != ENOENT)
        die("unlink '%s': %s", lock_path, strerror(errno));
    free(lock_path);
}

extern int cmd_sort_main(int argc, char **argv);  /* for --write-index */

static void maybe_report_import_progress(view_t *v,
                                         int64_t n_in,
                                         int64_t n_kept,
                                         int64_t n_missing,
                                         int64_t voffset,
                                         int64_t file_size_bytes) {
    lmdb_backend_commit(v->lmdb);
    lmdb_backend_begin(v->lmdb);
    if (file_size_bytes > 0 && voffset >= 0) {
        int64_t compressed_offset = voffset >> 16;
        if (compressed_offset < 0)
            compressed_offset = 0;
        if (compressed_offset > file_size_bytes)
            compressed_offset = file_size_bytes;
        int pct = (int)((compressed_offset * 100) / file_size_bytes);
        if (pct > 99)
            pct = 99;
        fprintf(stderr,
                "isamtools import: %d%%, processed %lld records (%lld kept, %lld missing)\r",
                pct,
                (long long)n_in,
                (long long)n_kept,
                (long long)n_missing);
        return;
    }
    fprintf(stderr,
            "isamtools import: processed %lld records (%lld kept, %lld missing)\r",
            (long long)n_in,
            (long long)n_kept,
            (long long)n_missing);
}

static int64_t import_batch_limit(int64_t n_in) {
    (void)n_in;
    return 10000;
}

int cmd_import(int argc, char **argv) {
    const char *store_path = NULL;
    const char *store_idx  = NULL;
    const char *aln_path   = NULL;
    const char *name_alias = NULL;
    const char *carry_tags = NULL;
    const char *keep_alignment_tags = NULL;
    const char *drop_tags = NULL;
    const char *drop_source_tags = NULL;
    int strict_tag_policy = 0;
    int verbose_tag_policy = 0;
    int verify_after_import = 0;
    int allow_missing = 0;
    int write_index   = 0;

    static struct option longopts[] = {
        {"store",         required_argument, 0, 1000},
        {"store-index",   required_argument, 0, 1001},
        {"aln",           required_argument, 0, 1002},
        {"name",          required_argument, 0, 1003},
        {"allow-missing", no_argument,       0, 1005},
        {"write-index",   no_argument,       0, 1006},
        {"carry-tags",    required_argument, 0, 1007},
        {"keep-alignment-tags", required_argument, 0, 1008},
        {"drop-tags",     required_argument, 0, 1009},
        {"drop-source-tags", required_argument, 0, 1010},
        {"strict-tag-policy", no_argument,   0, 1011},
        {"verbose-tag-policy", no_argument,  0, 1012},
        {"check",        no_argument,       0, 1013},
        {"help",          no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: store_path = optarg; break;
        case 1001: store_idx  = optarg; break;
        case 1002: aln_path   = optarg; break;
        case 1003: name_alias = optarg; break;
        case 1007: carry_tags = optarg; break;
        case 1005: allow_missing = 1;   break;
        case 1006: write_index   = 1;   break;
        case 1008: keep_alignment_tags = optarg; break;
        case 1009: drop_tags = optarg; break;
        case 1010: drop_source_tags = optarg; break;
        case 1011: strict_tag_policy = 1; break;
        case 1012: verbose_tag_policy = 1; break;
        case 1013: verify_after_import = 1; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    if (!store_path) die("--store is required");
    const char *input_path = argv[optind];

    char *aln_path_alloc = NULL;
    if (!aln_path) {
        aln_path_alloc = default_aln_path_for(input_path);
        if (!aln_path_alloc)
            die("--aln is required when importing from stdin");
        aln_path = aln_path_alloc;
    }
    char *final_aln_path = xstrdup(aln_path);

    char *store_idx_alloc = NULL;
    if (!store_idx) { store_idx_alloc = isam_index_path_for(store_path);
                      store_idx = store_idx_alloc; }

    /* Open store + verify QNAME index. */
    store_t *store = store_open(store_path, store_idx);
    if (!store_has_name_index(store))
        die("store index '%s' has no QNAME index; "
            "run `isamtools index %s` first", store_idx, store_path);

    /* Open input stream. */
    samFile *in = sam_open(input_path, "r");
    if (!in) die("cannot open input '%s'", input_path);
    const htsFormat *input_fmt = hts_get_format(in);
    int input_is_bam = (input_fmt && input_fmt->format == bam);
    sam_hdr_t *in_hdr = sam_hdr_read(in);
    if (!in_hdr) die("cannot read header from input '%s'", input_path);

    int64_t input_size_bytes = -1;
    if (input_is_bam && strcmp(input_path, "-") != 0) {
        struct stat input_st;
        if (stat(input_path, &input_st) == 0)
            input_size_bytes = (int64_t)input_st.st_size;
    }

    /* Open / create output view. */
    prepare_import_output_path(aln_path);
    view_t *v = view_open_write(aln_path);

    /* Persist view metadata. */
    view_meta_set(v, "backend",           "lmdb");
    view_meta_set(v, "isam_type",         "alignment");
    view_meta_set(v, "isam_version",      ISAM_SCHEMA_VIEW);
    view_meta_set(v, "isamtools_version", ISAM_VERSION);
    view_meta_set(v, "store_path",        store_path);
    view_meta_set(v, "store_index_path",  store_idx);
    char *store_name = lmdb_backend_meta_get(store->lmdb, "name");
    char *shh = lmdb_backend_meta_get(store->lmdb, "store_header_hash");
    char *ssh = lmdb_backend_meta_get(store->lmdb, "store_sq_hash");
    char *sml = lmdb_backend_meta_get(store->lmdb, "store_max_read_len");
    if (store_name) { view_meta_set(v, "store_name", store_name); free(store_name); }
    if (shh) { view_meta_set(v, "store_header_hash", shh); free(shh); }
    if (ssh) { view_meta_set(v, "store_sq_hash", ssh); free(ssh); }
    if (sml) { view_meta_set(v, "store_max_read_len", sml); free(sml); }
    if (name_alias) view_meta_set(v, "name", name_alias);
    char *store_shared_tags = lmdb_backend_meta_get(store->lmdb, "shared_tags");
    char *store_shared_tags_for_message =
        xstrdup(store_shared_tags ? store_shared_tags : "");
    if (store_shared_tags) {
        view_meta_set(v, "shared_tags_snapshot", store_shared_tags);
    }
    char *ts = now_iso8601(); view_meta_set(v, "created_at", ts); free(ts);
    view_meta_set(v, "indexed_coord",    "0");
    view_meta_set(v, "indexed_name",     "0");
    view_meta_set(v, "indexed_interval", "0");
    view_meta_set(v, "has_import_table", "1");
    view_meta_set(v, "import_state",     "in_progress");

        fprintf(stderr,
            "isamtools import: importing alignments from %s -> %s using store %s\n",
            input_path, aln_path, store_path);
        fprintf(stderr,
            "isamtools import: matching QNAMEs against %s%s%s\n",
            store_idx,
            allow_missing ? ", skipping missing reads" : ", requiring all reads in the store",
            write_index ? ", then building all secondary indexes" : "");
        fprintf(stderr,
            "isamtools import: aux tags: canonical-backed tags are restored from the store BAM/CRAM when missing; import differences and carried tags are stored once per read; alignment-local tags remain on each alignment record\n");

    view_set_header(v, in_hdr);

    isam_tag_classifier_t classifier;
    isam_tag_classifier_init(&classifier, in_hdr, store_shared_tags, NULL, NULL,
                             carry_tags, keep_alignment_tags);
    if (store_shared_tags) free(store_shared_tags);

    isam_tag_set_t dropped_tags;
    isam_tag_set_t source_dropped_tags;
    isam_tag_set_init(&dropped_tags);
    isam_tag_set_init(&source_dropped_tags);
    if (isam_tag_set_parse_csv(&dropped_tags, drop_tags) != 0)
        die("invalid drop tag list '%s'", drop_tags);
    if (isam_tag_set_parse_csv(&source_dropped_tags, drop_source_tags) != 0)
        die("invalid import-drop tag list '%s'", drop_source_tags);
    for (size_t i = 0; i < dropped_tags.n; ++i) {
        for (size_t j = 0; j < source_dropped_tags.n; ++j) {
            if (source_dropped_tags.items[j] != dropped_tags.items[i])
                continue;
            char tag[3];
            isam_tag_decode(dropped_tags.items[i], tag);
            warnx_("import: --drop-tags wins over --drop-source-tags for '%s'", tag);
            memmove(&source_dropped_tags.items[j],
                    &source_dropped_tags.items[j + 1],
                    (source_dropped_tags.n - j - 1) *
                        sizeof source_dropped_tags.items[0]);
            --source_dropped_tags.n;
            break;
        }
    }

    char *shared_csv = isam_tag_set_join_csv(&classifier.effective_shared_tags);
    char *carried_csv = isam_tag_set_join_csv(&classifier.effective_carried_tags);
    char *alignment_csv = isam_tag_set_join_csv(&classifier.forced_alignment_tags);
    char *dropped_csv = isam_tag_set_join_csv(&dropped_tags);
    char *source_dropped_csv = isam_tag_set_join_csv(&source_dropped_tags);
    char *carry_user_csv = isam_tag_set_join_csv(&classifier.configured_carried_tags);
    char *keep_alignment_user_csv = isam_tag_set_join_csv(&classifier.forced_alignment_tags);
    view_meta_set(v, "shared_tags_snapshot", shared_csv);
    view_meta_set(v, "carried_tags", carried_csv);
    view_meta_set(v, "alignment_tags", alignment_csv);
    view_meta_set(v, "dropped_tags", dropped_csv);
    view_meta_set(v, "source_dropped_tags", source_dropped_csv);
    view_meta_set(v, "carry_tags_user_added", carry_user_csv);
    view_meta_set(v, "keep_alignment_tags_user_added", keep_alignment_user_csv);
    view_meta_set(v, "tag_policy_version", "1");

    fprintf(stderr, "isamtools import: tag policy before scanning\n");
    fprintf(stderr, "  store file: %s\n", store_path);
    fprintf(stderr, "  import file: %s\n", input_path);
    fprintf(stderr, "  canonical-backed tags: %s\n",
            store_shared_tags_for_message[0] ? store_shared_tags_for_message : "none");
    fprintf(stderr, "    if absent from the import record, they are restored from the canonical BAM/CRAM during view\n");
    fprintf(stderr, "    if the import record has a different value, that override is stored once per read in read_aux_by_id\n");
    fprintf(stderr, "  extra import read-level tags stored once per read: %s\n",
            carried_csv[0] ? carried_csv : "none");
    fprintf(stderr, "  tags kept per alignment record: built-in alignment tags%s%s\n",
            keep_alignment_user_csv[0] ? ", plus " : "",
            keep_alignment_user_csv[0] ? keep_alignment_user_csv : "");
    fprintf(stderr, "  tags removed from final output (--drop-tags): %s\n",
            dropped_csv[0] ? dropped_csv : "none");
    fprintf(stderr, "  import tag copies ignored; canonical fallback still allowed (--drop-source-tags): %s\n",
            source_dropped_csv[0] ? source_dropped_csv : "none");
    fprintf(stderr, "  observed tags and storage counts will be printed after the scan\n");

    lmdb_backend_begin(v->lmdb);
    bam1_t *b = bam_init1();
    bam1_t *canon_b = bam_init1();
    payload_cache_t *canonical_cache = payload_cache_create(64 * 1024 * 1024);
    int64_t n_in = 0, n_kept = 0, n_missing = 0;
    int64_t max_ref_span = 0;
    int64_t max_query_len = 0;
    int64_t n_source_tag_fields_stored = 0;
    int64_t n_alignment_tag_fields_kept_inline = 0;
    int64_t n_dropped_tag_fields = 0;
    int64_t n_source_dropped_tag_fields = 0;
    int64_t n_unknown_tag_fields_kept_inline = 0;
    int64_t batch = 0;
    khash_t(u64set) *conflict_reads = kh_init(u64set);
    if (!conflict_reads)
        die("out of memory");
    khash_t(u64set) *seen_read_tags = kh_init(u64set);
    if (!seen_read_tags)
        die("out of memory");
    read_aux_dict_trainer_t read_aux_trainer = {0};
    tag_count_t *conflict_counts = NULL;
    size_t n_conflict_counts = 0, cap_conflict_counts = 0;
    isam_tag_set_t observed_alignment_tags;
    isam_tag_set_init(&observed_alignment_tags);
    isam_tag_set_t observed_unknown_tags;
    isam_tag_set_init(&observed_unknown_tags);
    isam_tag_set_t observed_source_overrides;
    isam_tag_set_init(&observed_source_overrides);
    isam_tag_set_t observed_carried_tags;
    isam_tag_set_init(&observed_carried_tags);
    isam_tag_set_t observed_drop_tags;
    isam_tag_set_init(&observed_drop_tags);
    isam_tag_set_t observed_source_drop_tags;
    isam_tag_set_init(&observed_source_drop_tags);
    while (1) {
        int64_t voffset = -1;
        if (input_is_bam)
            voffset = bgzf_tell(in->fp.bgzf);

        int rc = sam_read1(in, in_hdr, b);
        if (rc == -1) break;
        if (rc < -1)  die("error reading input '%s' at record %lld",
                          input_path, (long long)n_in);

        const char *qname = bam_get_qname(b);
        int64_t qh = qname_hash(qname);

        /* Lookup in canonical store. */
        store_record_t srec;
        int hit = store_lookup_qname(store, qname, &srec);
        if (!hit) {
            if (!allow_missing)
                die("QNAME '%s' not found in canonical store index '%s' "
                    "(canonical indexing prefers primary alignments and falls back to non-hard-clipped supplementary alignments; use --allow-missing to skip)",
                    qname, store_path);
            ++n_missing;
            ++n_in;
            continue;
        }

        int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
        int n_cigar = b->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(b);
        bam1_t *filtered = bam_dup1(b);
        if (!filtered) die("out of memory");
        int canonical_loaded = 0;
        uint8_t *read_aux = NULL;
        size_t read_aux_len = 0;
        size_t read_aux_cap = 0;
        int read_aux_changed = 0;
        if (lmdb_read_aux_get(v->lmdb, srec.read_id, &read_aux, &read_aux_len))
            read_aux_cap = read_aux_len;

        errno = 0;
        for (uint8_t *field = bam_aux_first(b); field; field = bam_aux_next(b, field)) {
            char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
            if (isam_tag_set_contains(&dropped_tags, tag)) {
                uint8_t *drop = bam_aux_get(filtered, tag);
                if (drop && bam_aux_del(filtered, drop) != 0)
                    die("failed to strip dropped aux tag '%s'", tag);
                isam_tag_set_add(&observed_drop_tags, tag);
                ++n_dropped_tag_fields;
                continue;
            }
            if (isam_tag_set_contains(&source_dropped_tags, tag)) {
                uint8_t *drop = bam_aux_get(filtered, tag);
                if (drop && bam_aux_del(filtered, drop) != 0)
                    die("failed to strip import-only dropped aux tag '%s'", tag);
                isam_tag_set_add(&observed_source_drop_tags, tag);
                ++n_source_dropped_tag_fields;
                continue;
            }

            isam_tag_scope_t scope = isam_tag_classifier_scope(&classifier, tag);
            if (scope == ISAM_TAG_SCOPE_SHARED || scope == ISAM_TAG_SCOPE_CARRIED) {
                int payload_len = isam_aux_field_payload_len(field);
                if (payload_len < 0)
                    die("invalid aux tag '%s' in '%s'", tag, input_path);
                char type[2] = { bam_aux_type(field), 0 };
                int store_has_same = 0;
                if (scope == ISAM_TAG_SCOPE_SHARED) {
                    if (!canonical_loaded) {
                        int frc = payload_cache_fetch(canonical_cache, store,
                                                      &srec, canon_b);
                        if (frc != 0)
                            die("payload fetch failed for '%s' (rc=%d)",
                                srec.qname, frc);
                        canonical_loaded = 1;
                    }
                    store_has_same = aux_field_value_matches(bam_aux_get(canon_b, tag),
                                                             type[0],
                                                             field + 1,
                                                             (size_t)payload_len);
                }
                int stored_import_value = 0;
                const uint8_t *existing = aux_blob_find_field(read_aux, read_aux_len, tag);
                int seen_store_value =
                    scope == ISAM_TAG_SCOPE_SHARED && !existing &&
                    read_tag_seen_contains(seen_read_tags, srec.read_id, tag);
                if (existing) {
                    int same = aux_field_value_matches(existing, type[0],
                                                       field + 1,
                                                       (size_t)payload_len);
                    if (!same &&
                        read_tag_seen_add(conflict_reads, srec.read_id, tag)) {
                        tag_count_inc(&conflict_counts, &n_conflict_counts,
                                      &cap_conflict_counts, tag);
                    }
                } else if (seen_store_value) {
                    if (!store_has_same &&
                        read_tag_seen_add(conflict_reads, srec.read_id, tag)) {
                        tag_count_inc(&conflict_counts, &n_conflict_counts,
                                      &cap_conflict_counts, tag);
                    }
                } else {
                    if (store_has_same && scope == ISAM_TAG_SCOPE_SHARED) {
                        read_tag_seen_add(seen_read_tags, srec.read_id, tag);
                    } else {
                        aux_blob_append_field(&read_aux, &read_aux_len,
                                              &read_aux_cap, field);
                        read_aux_changed = 1;
                        stored_import_value = 1;
                        ++n_source_tag_fields_stored;
                    }
                }
                if (scope == ISAM_TAG_SCOPE_SHARED && stored_import_value)
                    isam_tag_set_add(&observed_source_overrides, tag);
                else if (scope == ISAM_TAG_SCOPE_CARRIED)
                    isam_tag_set_add(&observed_carried_tags, tag);
                uint8_t *drop = bam_aux_get(filtered, tag);
                if (drop && bam_aux_del(filtered, drop) != 0)
                    die("failed to strip stored-once aux tag '%s'", tag);
            } else if (scope == ISAM_TAG_SCOPE_ALIGNMENT) {
                isam_tag_set_add(&observed_alignment_tags, tag);
                ++n_alignment_tag_fields_kept_inline;
            } else if (scope == ISAM_TAG_SCOPE_UNKNOWN) {
                if (isam_tag_set_add(&observed_alignment_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, input_path);
                if (isam_tag_set_add(&observed_unknown_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, input_path);
                ++n_unknown_tag_fields_kept_inline;
            }
        }
        if (errno && errno != ENOENT)
            die("corrupt aux data in '%s' at record %lld",
                input_path, (long long)n_in);

        const uint8_t *aux = bam_get_aux(filtered);
        int aux_len = bam_get_l_aux(filtered);
        if ((int64_t)b->core.l_qseq > max_query_len)
            max_query_len = b->core.l_qseq;
        if (!is_unmapped && b->core.tid >= 0 && b->core.pos >= 0) {
            int64_t span = (int64_t)bam_endpos(b) - (int64_t)b->core.pos;
            if (span > max_ref_span) max_ref_span = span;
        }

        lmdb_view_put_read(v->lmdb, srec.read_id, qname, qh);
        lmdb_view_record_t rec = {
            .aln_id = n_kept + 1,
            .store_read_id = srec.read_id,
            .flag = b->core.flag,
            .tid = b->core.tid,
            .pos0 = (is_unmapped || b->core.pos < 0) ? -1 : (int64_t)b->core.pos,
            .end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b),
            .mapq = b->core.qual,
            .n_cigar = n_cigar,
            .cigar = (uint32_t *)cigar,
            .mate_tid = b->core.mtid,
            .mate_pos0 = b->core.mpos < 0 ? -1 : (int64_t)b->core.mpos,
            .tlen = (int64_t)b->core.isize,
            .aux = (uint8_t *)aux,
            .aux_len = aux_len,
            .record_order = n_kept,
            .is_unmapped = is_unmapped
        };
        lmdb_view_put_import_record(v->lmdb, &rec);
        if (read_aux_changed) {
            lmdb_read_aux_put(v->lmdb, srec.read_id, read_aux, read_aux_len);
            maybe_train_read_aux_dict(v->lmdb, &read_aux_trainer,
                                      read_aux, read_aux_len);
        }
        free(read_aux);
        bam_destroy1(filtered);
        store_record_clear(&srec);

        ++n_kept;
        ++n_in;
        ++batch;
        if (batch >= import_batch_limit(n_in)) {
            batch = 0;
            maybe_report_import_progress(v, n_in, n_kept, n_missing,
                                         voffset, input_size_bytes);
        }
    }
    lmdb_backend_commit(v->lmdb);
    bam_destroy1(b);
    bam_destroy1(canon_b);
    payload_cache_destroy(canonical_cache);
    sam_hdr_destroy(in_hdr);
    sam_close(in);

    view_meta_set_i64(v, "view_max_ref_span", max_ref_span);
    view_meta_set_i64(v, "view_max_query_len", max_query_len);
    char *observed_alignment_csv = isam_tag_set_join_csv(&observed_alignment_tags);
    char *observed_unknown_csv = isam_tag_set_join_csv(&observed_unknown_tags);
    char *source_override_csv = isam_tag_set_join_csv(&observed_source_overrides);
    char *observed_carried_csv = isam_tag_set_join_csv(&observed_carried_tags);
    char *observed_drop_csv = isam_tag_set_join_csv(&observed_drop_tags);
    char *observed_source_drop_csv = isam_tag_set_join_csv(&observed_source_drop_tags);
    view_meta_set(v, "alignment_tags", observed_alignment_csv);
    view_meta_set(v, "source_override_tags", source_override_csv);
    view_meta_set(v, "dropped_tags_seen", observed_drop_csv);
    view_meta_set(v, "source_dropped_tags_seen", observed_source_drop_csv);
    view_meta_set_i64(v, "carried_tag_conflicts_seen", (int64_t)n_conflict_counts);
    if (n_unknown_tag_fields_kept_inline > 0)
        warnx_("import: unknown aux tags default to alignment retention: %s",
               observed_unknown_csv);
    if (observed_drop_csv[0])
        warnx_("import: dropped import fields for tags: %s", observed_drop_csv);
    if (observed_source_drop_csv[0])
        warnx_("import: ignored import-only fields for tags: %s", observed_source_drop_csv);
    if (n_conflict_counts > 0) {
        warnx_("import: stored-once tag conflicts found; first observed values were retained");
        for (size_t i = 0; i < n_conflict_counts; ++i) {
            char tag[3];
            isam_tag_decode(conflict_counts[i].tag_code, tag);
            warnx_("import: %s: %lld reads with conflicts",
                   tag, (long long)conflict_counts[i].count);
        }
        if (strict_tag_policy)
            die("import: stored-once tag conflicts are errors under --strict-tag-policy");
    }
    view_meta_set(v, "import_state", "complete");
    fprintf(stderr, verbose_tag_policy
            ? "isamtools import: observed aux tag summary\n"
            : "isamtools import: observed aux tags\n");
    fprintf(stderr, "  canonical-backed tags with import overrides stored once per read: %s\n",
            source_override_csv[0] ? source_override_csv : "none");
    fprintf(stderr, "  extra import tags stored once per read: %s\n",
            observed_carried_csv[0] ? observed_carried_csv : "none");
    fprintf(stderr, "  tags stored per alignment record: %s\n",
            observed_alignment_csv[0] ? observed_alignment_csv : "none");
    fprintf(stderr, "  unknown tags stored per alignment record: %s\n",
            observed_unknown_csv[0] ? observed_unknown_csv : "none");
    fprintf(stderr, "  tags removed from final output: %s\n",
            observed_drop_csv[0] ? observed_drop_csv : "none");
    fprintf(stderr, "  import tag copies ignored while allowing canonical fallback: %s\n",
            observed_source_drop_csv[0] ? observed_source_drop_csv : "none");
    fprintf(stderr, "  final-output rule: --drop-tags wins; import values win over canonical fallback; canonical fallback fills missing tags\n");
    free(observed_alignment_csv);
    free(observed_unknown_csv);
    free(source_override_csv);
    free(observed_carried_csv);
    free(observed_drop_csv);
    free(observed_source_drop_csv);

    view_close(v);
    store_close(store);
    free(store_idx_alloc);
    free(aln_path_alloc);
    isam_tag_set_destroy(&dropped_tags);
    isam_tag_set_destroy(&source_dropped_tags);
    isam_tag_set_destroy(&observed_alignment_tags);
    isam_tag_set_destroy(&observed_unknown_tags);
    isam_tag_set_destroy(&observed_source_overrides);
    isam_tag_set_destroy(&observed_carried_tags);
    isam_tag_set_destroy(&observed_drop_tags);
    isam_tag_set_destroy(&observed_source_drop_tags);
    isam_tag_classifier_destroy(&classifier);
    kh_destroy(u64set, conflict_reads);
    kh_destroy(u64set, seen_read_tags);
    free(read_aux_trainer.samples.s);
    free(conflict_counts);
    free(shared_csv);
    free(carried_csv);
    free(alignment_csv);
    free(dropped_csv);
    free(source_dropped_csv);
    free(carry_user_csv);
    free(keep_alignment_user_csv);
    free(store_shared_tags_for_message);

    fprintf(stderr,
            "isamtools import: completed %s\n",
            final_aln_path);
    fprintf(stderr,
            "isamtools import: processed %lld records (%lld kept, %lld missing, max query length %lld, max reference span %lld)\n",
            (long long)n_in,
            (long long)n_kept,
            (long long)n_missing,
            (long long)max_query_len,
            (long long)max_ref_span);
    fprintf(stderr, "isamtools import: aux tag storage summary\n");
    fprintf(stderr,
            "  read_aux_by_id fields written once per read: %lld (canonical overrides plus extra import read tags)\n",
            (long long)n_source_tag_fields_stored);
    fprintf(stderr,
            "  aln_by_id inline aux fields written per alignment: %lld (alignment-local tags; includes %lld unknown fields)\n",
            (long long)(n_alignment_tag_fields_kept_inline +
                        n_unknown_tag_fields_kept_inline),
            (long long)n_unknown_tag_fields_kept_inline);
    fprintf(stderr,
            "  import aux fields removed by --drop-tags: %lld\n",
            (long long)n_dropped_tag_fields);
    fprintf(stderr,
            "  import aux fields ignored by --drop-source-tags: %lld\n",
            (long long)n_source_dropped_tag_fields);

    if (verify_after_import)
        isam_check_view_against_bam(input_path, final_aln_path, NULL, NULL);

    if (write_index) {
        fprintf(stderr,
                "isamtools import: running post-import sort --all on %s\n",
                final_aln_path);
        char *args[] = { "sort", "--all", (char *)final_aln_path, NULL };
        int rc = cmd_sort_main(3, args);
        free(final_aln_path);
        return rc;
    }
    free(final_aln_path);
    return 0;
}
