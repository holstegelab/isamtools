#include "common.h"
#include "backend_lmdb.h"
#include "aux_tags.h"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/cram.h>
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
    "Build the canonical store index <store-root>.ssam.\n"
        "By default builds the QNAME index (required for `import`).\n"
        "\n"
        "Options:\n"
        "  --name NAME      Logical alias for this store (e.g. chm13).\n"
        "  --share-tags TAGS     Aux tags reconstructed from canonical BAM/CRAM when missing.\n"
        "  --no-share-tags TAGS  Aux tags excluded from canonical reconstruction.\n"
        "  --strict-tag-policy   Fail when unknown canonical tags are unresolved.\n"
        "  --verbose-tag-policy  Print a detailed tag policy summary.\n"
        "  --yes                 Accept tag policy warnings.\n"
    "  -o, --output P   Override default output path (<store-root>.ssam).\n"
        "  -T, --reference P  Reference FASTA (required for a CRAM canonical store).\n"
        "  -h, --help       Show this help.\n",
        stderr);
}

/* Scan a CRAM canonical store's container index and return, in container order,
 * each container's byte offset and record count. Requires a .crai. The view
 * side uses this table to fetch a read's payload by (container offset, index). */
static int64_t build_cram_container_table(const char *store_path,
                                          int64_t **coff_out, int32_t **nrec_out) {
    samFile *cf = sam_open(store_path, "r");
    if (!cf) die("cannot open '%s' for CRAM container scan", store_path);
    sam_hdr_t *ch = sam_hdr_read(cf);
    if (!ch) die("cannot read header of '%s'", store_path);
    hts_idx_t *idx = sam_index_load(cf, store_path);
    if (!idx)
        die("CRAM canonical store '%s' has no .crai index; run `samtools index` on it first",
            store_path);
    cram_fd *fd = cf->fp.cram;
    int64_t nc = cram_num_containers(fd);
    if (nc <= 0)
        die("cram_num_containers returned %lld for '%s'", (long long)nc, store_path);
    int64_t *coff = xmalloc((size_t)nc * sizeof *coff);
    int32_t *nrec = xmalloc((size_t)nc * sizeof *nrec);
    for (int64_t k = 0; k < nc; ++k) {
        off_t off = cram_container_num2offset(fd, k);
        if (off < 0)
            die("cram_container_num2offset(%lld) failed for '%s'", (long long)k, store_path);
        if (cram_seek(fd, off, SEEK_SET) < 0)
            die("cram_seek during container scan of '%s'", store_path);
        cram_container *c = cram_read_container(fd);
        if (!c)
            die("cram_read_container %lld failed for '%s'", (long long)k, store_path);
        coff[k] = (int64_t)off;
        nrec[k] = cram_container_get_num_records(c);
    }
    hts_idx_destroy(idx);
    sam_hdr_destroy(ch);
    sam_close(cf);
    *coff_out = coff;
    *nrec_out = nrec;
    return nc;
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
    int is_unmapped;
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
    int64_t n_indexed_mapped_reads;
    int64_t n_indexed_unmapped_reads;
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
                                            int is_supplementary,
                                            int is_unmapped) {
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
    kh_value(map, iter).is_unmapped = is_unmapped;
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

static uint64_t dbi_page_bytes(MDB_txn *txn, MDB_dbi dbi, size_t page_size) {
    if (!dbi)
        return 0;
    MDB_stat st;
    int rc = mdb_stat(txn, dbi, &st);
    if (rc)
        die("mdb_stat: %s", mdb_strerror(rc));
    uint64_t pages = (uint64_t)st.ms_branch_pages +
                     (uint64_t)st.ms_leaf_pages +
                     (uint64_t)st.ms_overflow_pages;
    return pages * (uint64_t)page_size;
}

static uint64_t lmdb_page_bytes(lmdb_backend_t *lmdb) {
    MDB_stat env_st;
    int rc = mdb_env_stat(lmdb->env, &env_st);
    if (rc)
        die("mdb_env_stat: %s", mdb_strerror(rc));

    MDB_txn *txn = NULL;
    rc = mdb_txn_begin(lmdb->env, NULL, MDB_RDONLY, &txn);
    if (rc)
        die("txn_begin: %s", mdb_strerror(rc));
    uint64_t bytes = 0;
    bytes += dbi_page_bytes(txn, lmdb->meta, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->canonical_by_id, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->canonical_by_qname, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->read_by_id, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->read_by_qname, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_by_id, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_import_by_order, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_by_coord, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_unmapped, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_by_read, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_by_bin, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->aln_checkpoint, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->subset_by_name, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->subset_read, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->read_aux_by_id, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->read_tag, env_st.ms_psize);
    bytes += dbi_page_bytes(txn, lmdb->annotation_meta, env_st.ms_psize);
    mdb_txn_abort(txn);
    return bytes;
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
    char *header_shared_tags = isam_tag_set_join_csv(&classifier.header_shared_tags);
    char *effective_shared_tags = isam_tag_set_join_csv(&classifier.effective_shared_tags);
    lmdb_store_put_meta(lmdb, "no_share_tags", no_shared_tags);
    lmdb_store_put_meta(lmdb, "share_tags_user_added", user_shared_tags);
    lmdb_store_put_meta(lmdb, "tag_policy_version", "1");
    fprintf(stderr, "isamtools index: tag policy before scanning\n");
    fprintf(stderr, "  store file: %s\n", store_path);
    fprintf(stderr, "  canonical aux storage: not copied into the .ssam; reconstructed from the canonical BAM/CRAM payload when needed\n");
    fprintf(stderr, "  read-level reconstruction candidates from policy/header: %s\n",
            effective_shared_tags[0] ? effective_shared_tags : "none");
    fprintf(stderr, "  header-declared read-level candidates: %s\n",
            header_shared_tags[0] ? header_shared_tags : "none");
    fprintf(stderr, "  tags explicitly added by user: %s\n",
            user_shared_tags[0] ? user_shared_tags : "none");
    fprintf(stderr, "  canonical alignment-only tags stay only in the source BAM/CRAM payload\n");
    fprintf(stderr, "  tags excluded from sharing by user request: %s\n",
            no_shared_tags[0] ? no_shared_tags : "none");
    fprintf(stderr, "  note: SAM/BAM headers do not generally enumerate every aux tag; observed tag statistics are printed after the scan\n");

    bam1_t *b = bam_init1();
    int64_t max_read_len = 0;
    int have_voffset = is_bam;

    /* CRAM canonical: build the container table and cumulative record bases so
     * each read can be tagged with a (container offset, within-container index)
     * payload locator, the CRAM analogue of the BAM BGZF virtual offset. */
    int64_t *cram_coff = NULL;
    int32_t *cram_nrec = NULL;
    int64_t *cram_base = NULL;
    int64_t  cram_nc = 0;
    int64_t  cram_ci = 0;
    if (!is_bam) {
        cram_nc = build_cram_container_table(store_path, &cram_coff, &cram_nrec);
        cram_base = xmalloc((size_t)(cram_nc + 1) * sizeof *cram_base);
        int64_t cum = 0;
        for (int64_t k = 0; k < cram_nc; ++k) { cram_base[k] = cum; cum += cram_nrec[k]; }
        cram_base[cram_nc] = cum;
        lmdb_store_put_cram_container_table(lmdb, cram_coff, cram_nrec, cram_nc);
        fprintf(stderr,
                "isamtools index: CRAM canonical: %lld containers, %lld records "
                "(payload located by container offset + within-container index)\n",
                (long long)cram_nc, (long long)cum);
    }

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

        /* Locate this raw record within the CRAM container table. The raw index
         * counts every record (incl. secondary/supplementary) so it stays in
         * lockstep with the container record counts and a sequential rescan. */
        int64_t cram_loc_coffset = -1;
        int32_t cram_loc_within = 0;
        int     cram_has_loc = 0;
        if (!is_bam) {
            int64_t raw_i = stats.n_alignment_records - 1;
            while (cram_ci + 1 < cram_nc && raw_i >= cram_base[cram_ci + 1]) cram_ci++;
            if (raw_i < cram_base[cram_ci] || raw_i >= cram_base[cram_ci + 1])
                die("cram record %lld outside container table [%lld,%lld)",
                    (long long)raw_i, (long long)cram_base[cram_ci],
                    (long long)cram_base[cram_ci + 1]);
            cram_loc_coffset = cram_coff[cram_ci];
            cram_loc_within = (int32_t)(raw_i - cram_base[cram_ci]);
            cram_has_loc = 1;
        }

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
                .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0,
                .cram_coffset = cram_loc_coffset,
                .cram_within = cram_loc_within,
                .has_cram_loc = cram_has_loc
            };
            lmdb_store_put_record(lmdb, &rec);
            canonical_map_put(canonical_reads, qname, read_id, 1, is_unmapped);
            ++stats.n_indexed_records;
            if (is_unmapped) ++stats.n_indexed_unmapped_reads;
            else ++stats.n_indexed_mapped_reads;
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

            int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
            if (state && state->is_supplementary) {
                ++stats.n_primary_upgrades;
                --stats.n_active_supplementary_records;
                if (state->is_unmapped) --stats.n_indexed_unmapped_reads;
                else --stats.n_indexed_mapped_reads;
                if (is_unmapped) ++stats.n_indexed_unmapped_reads;
                else ++stats.n_indexed_mapped_reads;
            }
            strset_remove(hardclip_only_qnames, qname);

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
                .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0,
                .cram_coffset = cram_loc_coffset,
                .cram_within = cram_loc_within,
                .has_cram_loc = cram_has_loc
            };
            if (state && state->is_supplementary) {
                lmdb_store_update_record(lmdb, &rec);
                state->is_supplementary = 0;
                state->is_unmapped = is_unmapped;
            } else {
                lmdb_store_put_record(lmdb, &rec);
                canonical_map_put(canonical_reads, qname, read_id, 0, is_unmapped);
                ++stats.n_indexed_records;
                if (is_unmapped) ++stats.n_indexed_unmapped_reads;
                else ++stats.n_indexed_mapped_reads;
            }
        }

        if ((int64_t)b->core.l_qseq > max_read_len)
            max_read_len = b->core.l_qseq;

        maybe_report_progress(lmdb, &stats, voffset, stbuf->st_size);
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(input);
    free(cram_coff);
    free(cram_nrec);
    free(cram_base);

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
    uint64_t page_bytes = lmdb_page_bytes(lmdb);

    stats.n_supplementary_fallback_reads = lmdb_store_count_supplementary_records(lmdb);
    if (stats.n_supplementary_fallback_reads != stats.n_active_supplementary_records)
        die("internal error: supplementary fallback count mismatch for '%s'",
            store_path);
    stats.n_supplementary_hard_clip_only_reads =
        (int64_t)kh_size(hardclip_only_qnames);
    lmdb_backend_close(lmdb);

    if (stats.n_alignment_records >= 10000)
        fputc('\n', stderr);
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
            ? "isamtools index: observed canonical aux statistics\n"
            : "isamtools index: observed canonical aux stats\n");
    fprintf(stderr, "  read-level tags observed and reconstructible from canonical BAM/CRAM: %s\n",
            shared_tags[0] ? shared_tags : "none");
    fprintf(stderr, "  alignment-only tags observed only in canonical BAM/CRAM payload: %s\n",
            observed_alignment_csv[0] ? observed_alignment_csv : "none");
    fprintf(stderr, "  tags actually excluded from sharing by user request: %s\n",
            excluded_seen_csv[0] ? excluded_seen_csv : "none");
    fprintf(stderr, "  unknown tags observed only in canonical BAM/CRAM payload: %s\n",
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
    free(header_shared_tags);
    free(effective_shared_tags);
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
    fprintf(stderr,
            "isamtools index: canonical reads by mapping state: %lld mapped, %lld unmapped\n",
            (long long)stats.n_indexed_mapped_reads,
            (long long)stats.n_indexed_unmapped_reads);
    fprintf(stderr,
            "isamtools index: LMDB storage page bytes: %" PRIu64 "\n",
            page_bytes);
    return 0;
}

int cmd_index(int argc, char **argv) {
    const char *name_alias = NULL;
    const char *out_path   = NULL;
    const char *share_tags = NULL;
    const char *no_share_tags = NULL;
    const char *reference = NULL;
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
        {"reference", required_argument, 0, 'T'},
        {"help",   no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:T:h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: name_alias = optarg; break;
        case 1002: share_tags = optarg; break;
        case 1003: no_share_tags = optarg; break;
        case 1004: yes = 1; break;
        case 1005: strict_tag_policy = 1; break;
        case 1006: verbose_tag_policy = 1; break;
        case 'o':  out_path   = optarg; break;
        case 'T':  reference  = optarg; break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    const char *store_path = argv[optind];

    char *out_alloc = NULL;
    if (!out_path) {
        out_alloc = isam_store_index_path_for(store_path);
        out_path = out_alloc;
    }

    samFile *input = sam_open(store_path, "r");
    if (!input) die("cannot open '%s'", store_path);
    const htsFormat *fmt = hts_get_format(input);
    int is_bam = (fmt && fmt->format == bam);
    int is_cram = (fmt && fmt->format == cram);
    if (!is_bam && !is_cram)
        die("unsupported canonical store format (need BAM or CRAM)");

    if (is_cram) {
        if (!reference)
            die("indexing a CRAM canonical store requires -T/--reference");
        if (hts_set_fai_filename(input, reference) != 0)
            die("cannot set reference FASTA '%s' for '%s'", reference, store_path);
    }

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
