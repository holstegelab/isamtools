#include "common.h"
#include "backend_lmdb.h"
#include "aux_tags.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>
#include <inttypes.h>

static void usage(void) {
    fputs(
        "Usage: isamtools index [options] <store.bam|store.cram>\n"
        "\n"
    "Build the canonical store index <store-root>.isam.\n"
        "By default builds the QNAME index (required for `import`).\n"
        "\n"
        "Options:\n"
        "  --name NAME      Logical alias for this store (e.g. chm13).\n"
        "  --share-tags TAGS     Comma-separated aux tags shared from the canonical store.\n"
        "  --no-share-tags TAGS  Comma-separated aux tags excluded from canonical sharing.\n"
        "  --strict-tag-policy   Fail when unknown canonical tags are unresolved.\n"
        "  --verbose-tag-policy  Print a detailed tag policy summary.\n"
        "  --yes                 Accept tag policy warnings.\n"
    "  -o, --output P   Override default output path (<store-root>.isam).\n"
        "  -h, --help       Show this help.\n",
        stderr);
}

static void lmdb_set_header_meta(lmdb_backend_t *lmdb, sam_hdr_t *hdr) {
    const char *htxt = sam_hdr_str(hdr);
    lmdb_store_put_meta_i64(lmdb, "store_header_hash",
                            header_text_hash(htxt ? htxt : ""));

    kstring_t sq = { 0, 0, NULL };
    int n = sam_hdr_nref(hdr);
    for (int i = 0; i < n; ++i) {
        char lbuf[32];
        snprintf(lbuf, sizeof lbuf, "%" PRId64,
                 (int64_t)sam_hdr_tid2len(hdr, i));
        kputs(sam_hdr_tid2name(hdr, i), &sq);
        kputc('\t', &sq);
        kputs(lbuf, &sq);
        kputc('\n', &sq);
    }
    lmdb_store_put_meta_i64(lmdb, "store_sq_hash",
                            header_text_hash(sq.s ? sq.s : ""));
    free(sq.s);
}

static int bam_has_hard_clip(const bam1_t *b) {
    int n_cigar = b->core.n_cigar;
    const uint32_t *cigar = bam_get_cigar(b);
    for (int i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)
            return 1;
    }
    return 0;
}

KHASH_SET_INIT_STR(strset)

typedef struct canonical_state {
    int64_t read_id;
    int is_supplementary;
} canonical_state_t;

KHASH_MAP_INIT_STR(canonical_map, canonical_state_t)

typedef struct index_stats {
    int64_t n_alignment_records;
    int64_t n_primary_records;
    int64_t n_secondary_records;
    int64_t n_supplementary_records;
    int64_t n_unmapped_records;
    int64_t n_supplementary_hard_clipped_records;
    int64_t n_supplementary_fallback_reads;
    int64_t n_supplementary_hard_clip_only_reads;
    int64_t n_indexed_records;
    int64_t n_primary_upgrades;
    int64_t n_active_supplementary_records;
} index_stats_t;

static khash_t(strset) *strset_create(void) {
    khash_t(strset) *set = kh_init(strset);
    if (!set)
        die("out of memory");
    return set;
}

static int strset_contains(const khash_t(strset) *set, const char *key) {
    return kh_get(strset, set, key) != kh_end(set);
}

static void strset_add(khash_t(strset) *set, const char *key) {
    if (strset_contains(set, key))
        return;
    char *owned = xstrdup(key);
    int absent = 0;
    khiter_t iter = kh_put(strset, set, owned, &absent);
    if (absent < 0) {
        free(owned);
        die("out of memory");
    }
    (void)iter;
}

static int strset_remove(khash_t(strset) *set, const char *key) {
    khiter_t iter = kh_get(strset, set, key);
    if (iter == kh_end(set))
        return 0;
    free((char *)kh_key(set, iter));
    kh_del(strset, set, iter);
    return 1;
}

static void strset_destroy(khash_t(strset) *set) {
    if (!set)
        return;
    for (khiter_t iter = kh_begin(set); iter != kh_end(set); ++iter) {
        if (!kh_exist(set, iter))
            continue;
        free((char *)kh_key(set, iter));
    }
    kh_destroy(strset, set);
}

static khash_t(canonical_map) *canonical_map_create(void) {
    khash_t(canonical_map) *map = kh_init(canonical_map);
    if (!map)
        die("out of memory");
    return map;
}

static canonical_state_t *canonical_map_get(khash_t(canonical_map) *map,
                                            const char *qname) {
    khiter_t iter = kh_get(canonical_map, map, qname);
    if (iter == kh_end(map))
        return NULL;
    return &kh_value(map, iter);
}

static canonical_state_t *canonical_map_put(khash_t(canonical_map) *map,
                                            const char *qname,
                                            int64_t read_id,
                                            int is_supplementary) {
    khiter_t iter = kh_get(canonical_map, map, qname);
    if (iter == kh_end(map)) {
        char *owned = xstrdup(qname);
        int absent = 0;
        iter = kh_put(canonical_map, map, owned, &absent);
        if (absent < 0) {
            free(owned);
            die("out of memory");
        }
    }
    kh_value(map, iter).read_id = read_id;
    kh_value(map, iter).is_supplementary = is_supplementary;
    return &kh_value(map, iter);
}

static void canonical_map_destroy(khash_t(canonical_map) *map) {
    if (!map)
        return;
    for (khiter_t iter = kh_begin(map); iter != kh_end(map); ++iter) {
        if (!kh_exist(map, iter))
            continue;
        free((char *)kh_key(map, iter));
    }
    kh_destroy(canonical_map, map);
}

static void maybe_report_progress(lmdb_backend_t *lmdb,
                                  index_stats_t *stats,
                                  int64_t voffset,
                                  int64_t file_size_bytes) {
    const int64_t first_report = 10000;
    const int64_t steady_report = 100000;

    if (stats->n_alignment_records == 0)
        return;
    if (stats->n_alignment_records < first_report ||
        (stats->n_alignment_records != first_report &&
         stats->n_alignment_records % steady_report != 0))
        return;
    lmdb_backend_commit(lmdb);
    lmdb_backend_begin(lmdb);
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
                "isamtools index: %d%%, processed %lld alignments, wrote %lld canonical reads (%lld supplementary fallback still active, %lld promoted to primary)\r",
                pct,
                (long long)stats->n_alignment_records,
                (long long)stats->n_indexed_records,
                (long long)stats->n_active_supplementary_records,
                (long long)stats->n_primary_upgrades);
        return;
    }
    fprintf(stderr,
            "isamtools index: processed %lld alignments, wrote %lld canonical reads (%lld supplementary fallback still active, %lld promoted to primary)\r",
            (long long)stats->n_alignment_records,
            (long long)stats->n_indexed_records,
            (long long)stats->n_active_supplementary_records,
            (long long)stats->n_primary_upgrades);
}

static int index_lmdb(const char *store_path, const char *out_path,
                      const char *name_alias, samFile *input,
                      sam_hdr_t *hdr, int is_bam,
                      const char *share_tags,
                      const char *no_share_tags,
                      int yes,
                      int strict_tag_policy,
                      int verbose_tag_policy,
                      const struct stat *stbuf) {
    struct stat out_st;
    if (stat(out_path, &out_st) == 0)
        die("LMDB output path '%s' already exists; choose a new path", out_path);
    if (errno != ENOENT)
        die("stat '%s': %s", out_path, strerror(errno));

    khash_t(canonical_map) *canonical_reads = canonical_map_create();
    khash_t(strset) *hardclip_only_qnames = strset_create();
    index_stats_t stats = {0};

    fprintf(stderr,
            "isamtools index: building canonical store index from %s -> %s\n",
            store_path, out_path);
    fprintf(stderr,
            "isamtools index: streaming canonical reads (prefer primary, fallback to non-hard-clipped supplementary)\n");

    lmdb_backend_t *lmdb = lmdb_backend_open(out_path, /*readonly*/0, 0);
    lmdb_backend_begin(lmdb);
    lmdb_store_put_meta(lmdb, "backend", "lmdb");
    lmdb_store_put_meta(lmdb, "isam_type", "store");
    lmdb_store_put_meta(lmdb, "isam_version", ISAM_SCHEMA_STORE);
    lmdb_store_put_meta(lmdb, "isamtools_version", ISAM_VERSION);
    lmdb_store_put_meta(lmdb, "store_path", store_path);
    lmdb_store_put_meta(lmdb, "store_format", is_bam ? "bam" : "cram");
    lmdb_store_put_meta_i64(lmdb, "store_size_bytes", (int64_t)stbuf->st_size);
    lmdb_store_put_meta_i64(lmdb, "store_mtime_ns",
                            (int64_t)stbuf->st_mtime * 1000000000LL);
    if (name_alias) lmdb_store_put_meta(lmdb, "name", name_alias);
    char *ts = now_iso8601();
    lmdb_store_put_meta(lmdb, "created_at", ts);
    free(ts);
    lmdb_set_header_meta(lmdb, hdr);

    isam_tag_classifier_t classifier;
    isam_tag_classifier_init(&classifier, hdr, NULL, share_tags, no_share_tags,
                             NULL, NULL);
    char *no_shared_tags = isam_tag_set_join_csv(&classifier.no_share_tags);
    char *user_shared_tags = isam_tag_set_join_csv(&classifier.configured_shared_tags);
    lmdb_store_put_meta(lmdb, "no_share_tags", no_shared_tags);
    lmdb_store_put_meta(lmdb, "share_tags_user_added", user_shared_tags);
    lmdb_store_put_meta(lmdb, "tag_policy_version", "1");
    fprintf(stderr, "isamtools index: tag policy before scanning\n");
    fprintf(stderr, "  store file: %s\n", store_path);
    fprintf(stderr, "  canonical read-level tags will be shared when observed\n");
    fprintf(stderr, "  canonical alignment-only tags stay in the source BAM/CRAM payload\n");
    fprintf(stderr, "  tags excluded from sharing by user request: %s\n",
            no_shared_tags[0] ? no_shared_tags : "none");
    fprintf(stderr, "  final observed tag summary will be printed after the scan\n");

    bam1_t *b = bam_init1();
    int64_t max_read_len = 0;
    int have_voffset = is_bam;
    isam_tag_set_t observed_shared_tags;
    isam_tag_set_t observed_alignment_tags;
    isam_tag_set_t unknown_canonical_tags;
    isam_tag_set_t excluded_seen_tags;
    isam_tag_set_init(&observed_shared_tags);
    isam_tag_set_init(&observed_alignment_tags);
    isam_tag_set_init(&unknown_canonical_tags);
    isam_tag_set_init(&excluded_seen_tags);

    while (1) {
        int64_t voffset = -1;
        if (have_voffset)
            voffset = bgzf_tell(input->fp.bgzf);

        int rc = sam_read1(input, hdr, b);
        if (rc == -1)
            break;
        if (rc < -1)
            die("error reading '%s' at record %lld",
                store_path, (long long)stats.n_alignment_records);

        const char *qname = bam_get_qname(b);
        ++stats.n_alignment_records;
        errno = 0;
        for (uint8_t *field = bam_aux_first(b); field; field = bam_aux_next(b, field)) {
            char tag[3] = { bam_aux_tag(field)[0], bam_aux_tag(field)[1], 0 };
            if (isam_tag_set_contains(&classifier.no_share_tags, tag)) {
                if (isam_tag_set_add(&excluded_seen_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, store_path);
                continue;
            }
            isam_tag_scope_t scope = isam_tag_classifier_scope(&classifier, tag);
            if (scope == ISAM_TAG_SCOPE_SHARED) {
                if (isam_tag_set_add(&observed_shared_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, store_path);
            } else if (scope == ISAM_TAG_SCOPE_ALIGNMENT) {
                if (isam_tag_set_add(&observed_alignment_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, store_path);
            } else if (scope == ISAM_TAG_SCOPE_UNKNOWN) {
                if (isam_tag_set_add(&unknown_canonical_tags, tag) != 0)
                    die("invalid aux tag '%s' in '%s'", tag, store_path);
            }
        }
        if (errno && errno != ENOENT)
            die("corrupt aux data in '%s' at record %lld",
                store_path, (long long)stats.n_alignment_records);

        if (b->core.flag & BAM_FUNMAP)
            ++stats.n_unmapped_records;

        if (b->core.flag & BAM_FSECONDARY) {
            ++stats.n_secondary_records;
            maybe_report_progress(lmdb, &stats, voffset, stbuf->st_size);
            continue;
        }

        canonical_state_t *state = canonical_map_get(canonical_reads, qname);
        if (b->core.flag & BAM_FSUPPLEMENTARY) {
            ++stats.n_supplementary_records;
            if (bam_has_hard_clip(b)) {
                ++stats.n_supplementary_hard_clipped_records;
                if (!state)
                    strset_add(hardclip_only_qnames, qname);
                maybe_report_progress(lmdb, &stats, voffset, stbuf->st_size);
                continue;
            }
            if (state) {
                maybe_report_progress(lmdb, &stats, voffset, stbuf->st_size);
                continue;
            }

            int64_t read_id = stats.n_indexed_records + 1;
            strset_remove(hardclip_only_qnames, qname);

            int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
            lmdb_store_record_t rec = {
                .read_id = read_id,
                .qname = (char *)qname,
                .qname_hash = qname_hash(qname),
                .tid = b->core.tid,
                .pos0 = (is_unmapped || b->core.pos < 0) ? -1 : (int64_t)b->core.pos,
                .end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b),
                .flag = b->core.flag,
                .is_unmapped = is_unmapped,
                .bam_voffset = (have_voffset && voffset >= 0) ? voffset : -1,
                .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0
            };
            lmdb_store_put_record(lmdb, &rec);
            canonical_map_put(canonical_reads, qname, read_id, 1);
            ++stats.n_indexed_records;
            ++stats.n_active_supplementary_records;
        } else {
            ++stats.n_primary_records;
            if (bam_has_hard_clip(b))
                die("primary canonical record '%s' in '%s' has hard clipping; "
                    "canonical store reads must keep full SEQ/QUAL",
                    qname, store_path);

            int64_t read_id = state ? state->read_id : (stats.n_indexed_records + 1);
            if (state && !state->is_supplementary)
                die("multiple primary canonical records for '%s' in '%s'",
                    qname, store_path);

            if (state && state->is_supplementary) {
                ++stats.n_primary_upgrades;
                --stats.n_active_supplementary_records;
            }
            strset_remove(hardclip_only_qnames, qname);

            int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
            lmdb_store_record_t rec = {
                .read_id = read_id,
                .qname = (char *)qname,
                .qname_hash = qname_hash(qname),
                .tid = b->core.tid,
                .pos0 = (is_unmapped || b->core.pos < 0) ? -1 : (int64_t)b->core.pos,
                .end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b),
                .flag = b->core.flag,
                .is_unmapped = is_unmapped,
                .bam_voffset = (have_voffset && voffset >= 0) ? voffset : -1,
                .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0
            };
            if (state && state->is_supplementary) {
                lmdb_store_update_record(lmdb, &rec);
                state->is_supplementary = 0;
            } else {
                lmdb_store_put_record(lmdb, &rec);
                canonical_map_put(canonical_reads, qname, read_id, 0);
                ++stats.n_indexed_records;
            }
        }

        if ((int64_t)b->core.l_qseq > max_read_len)
            max_read_len = b->core.l_qseq;

        maybe_report_progress(lmdb, &stats, voffset, stbuf->st_size);
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(input);

    lmdb_store_put_meta_i64(lmdb, "indexed_name", 1);
    lmdb_store_put_meta_i64(lmdb, "indexed_coord", 1);
    lmdb_store_put_meta_i64(lmdb, "store_max_read_len", max_read_len);
    char *shared_tags = isam_tag_set_join_csv(&observed_shared_tags);
    char *observed_alignment_csv = isam_tag_set_join_csv(&observed_alignment_tags);
    char *unknown_csv = isam_tag_set_join_csv(&unknown_canonical_tags);
    char *excluded_seen_csv = isam_tag_set_join_csv(&excluded_seen_tags);
    lmdb_store_put_meta(lmdb, "shared_tags", shared_tags);
    lmdb_store_put_meta(lmdb, "unknown_canonical_tags_seen", unknown_csv);
    if (unknown_csv[0] && strict_tag_policy && !yes)
        die("index: unknown canonical aux tags are not shared: %s", unknown_csv);
    lmdb_backend_commit(lmdb);

    stats.n_supplementary_fallback_reads = lmdb_store_count_supplementary_records(lmdb);
    if (stats.n_supplementary_fallback_reads != stats.n_active_supplementary_records)
        die("internal error: supplementary fallback count mismatch for '%s'",
            store_path);
    stats.n_supplementary_hard_clip_only_reads =
        (int64_t)kh_size(hardclip_only_qnames);
    lmdb_backend_close(lmdb);

    if (stats.n_supplementary_fallback_reads > 0)
        warnx_("index: used supplementary alignments as canonical records for %lld reads without a primary alignment",
               (long long)stats.n_supplementary_fallback_reads);
    if (stats.n_supplementary_hard_clip_only_reads > 0)
        warnx_("index: %lld reads had only hard-clipped supplementary alignments and were not indexed",
               (long long)stats.n_supplementary_hard_clip_only_reads);
    if (unknown_csv[0]) {
        warnx_("index: unknown canonical aux tags are not shared: %s", unknown_csv);
    }
    if (excluded_seen_csv[0])
        warnx_("index: canonical aux tags excluded from sharing by user request: %s",
               excluded_seen_csv);
    fprintf(stderr, verbose_tag_policy
            ? "isamtools index: observed tag policy summary\n"
            : "isamtools index: observed tag policy\n");
    fprintf(stderr, "  canonical read-level tags shared from store: %s\n",
            shared_tags[0] ? shared_tags : "none");
    fprintf(stderr, "  canonical alignment-only tags left in BAM/CRAM payload: %s\n",
            observed_alignment_csv[0] ? observed_alignment_csv : "none");
    fprintf(stderr, "  tags actually excluded from sharing by user request: %s\n",
            excluded_seen_csv[0] ? excluded_seen_csv : "none");
    fprintf(stderr, "  unknown tags left in BAM/CRAM payload: %s\n",
            unknown_csv[0] ? unknown_csv : "none");
    if (unknown_csv[0])
        fprintf(stderr, "  hint: use --share-tags %s to share unknown canonical tags\n",
                unknown_csv);
    isam_tag_classifier_destroy(&classifier);
    isam_tag_set_destroy(&observed_shared_tags);
    isam_tag_set_destroy(&observed_alignment_tags);
    isam_tag_set_destroy(&unknown_canonical_tags);
    isam_tag_set_destroy(&excluded_seen_tags);
    canonical_map_destroy(canonical_reads);
    strset_destroy(hardclip_only_qnames);
    free(shared_tags);
    free(observed_alignment_csv);
    free(no_shared_tags);
    free(user_shared_tags);
    free(unknown_csv);
    free(excluded_seen_csv);

    fprintf(stderr,
            "isamtools index: completed %s\n",
            out_path);
    fprintf(stderr,
            "isamtools index: scanned %lld alignments (%lld primary, %lld supplementary, %lld secondary, %lld unmapped)\n",
            (long long)stats.n_alignment_records,
            (long long)stats.n_primary_records,
            (long long)stats.n_supplementary_records,
            (long long)stats.n_secondary_records,
            (long long)stats.n_unmapped_records);
    fprintf(stderr,
            "isamtools index: indexed %lld canonical reads (%lld primary, %lld supplementary fallback, %lld supplementary promoted to primary, max read length %lld)\n",
            (long long)stats.n_indexed_records,
            (long long)(stats.n_indexed_records - stats.n_supplementary_fallback_reads),
            (long long)stats.n_supplementary_fallback_reads,
            (long long)stats.n_primary_upgrades,
            (long long)max_read_len);
    return 0;
}

int cmd_index(int argc, char **argv) {
    const char *name_alias = NULL;
    const char *out_path   = NULL;
    const char *share_tags = NULL;
    const char *no_share_tags = NULL;
    int yes = 0;
    int strict_tag_policy = 0;
    int verbose_tag_policy = 0;

    static struct option longopts[] = {
        {"name",   required_argument, 0, 1000},
        {"share-tags", required_argument, 0, 1002},
        {"no-share-tags", required_argument, 0, 1003},
        {"yes", no_argument, 0, 1004},
        {"strict-tag-policy", no_argument, 0, 1005},
        {"verbose-tag-policy", no_argument, 0, 1006},
        {"output", required_argument, 0, 'o'},
        {"help",   no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: name_alias = optarg; break;
        case 1002: share_tags = optarg; break;
        case 1003: no_share_tags = optarg; break;
        case 1004: yes = 1; break;
        case 1005: strict_tag_policy = 1; break;
        case 1006: verbose_tag_policy = 1; break;
        case 'o':  out_path   = optarg; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *store_path = argv[optind];

    char *out_alloc = NULL;
    if (!out_path) {
        out_alloc = isam_index_path_for(store_path);
        out_path = out_alloc;
    }

    samFile *input = sam_open(store_path, "r");
    if (!input) die("cannot open '%s'", store_path);
    const htsFormat *fmt = hts_get_format(input);
    int is_bam = (fmt && fmt->format == bam);
    int is_cram = (fmt && fmt->format == cram);
    if (!is_bam && !is_cram)
        die("unsupported canonical store format (need BAM or CRAM)");

    sam_hdr_t *hdr = sam_hdr_read(input);
    if (!hdr) die("cannot read header of '%s'", store_path);

    struct stat stbuf;
    if (stat(store_path, &stbuf) != 0)
        die("stat '%s': %s", store_path, strerror(errno));

    int rc = index_lmdb(store_path, out_path, name_alias, input, hdr,
                        is_bam, share_tags, no_share_tags, yes,
                        strict_tag_policy, verbose_tag_policy, &stbuf);
    free(out_alloc);
    return rc;
}
