#include "common.h"
#include "backend_lmdb.h"
#include "aux_tags.h"

#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <zstd.h>

typedef struct tag_stat {
    uint16_t tag_code;
    int64_t entries;
    int64_t payload_bytes;
    int64_t total_bytes;
} tag_stat_t;

static void usage(void) {
    fputs(
        "Usage: isamtools store_stats [options] <file.isam|file.ssam|store.bam|store.cram>\n"
        "\n"
        "Report LMDB table storage and logical record contents.\n"
        "\n"
        "Options:\n"
        "  --top N         Show only the top N tags per table (default all).\n"
        "  --storage-only  Skip logical scans; report only LMDB page usage.\n"
        "  -h, --help      Show this help.\n",
        stderr);
}

static char *resolve_isam_path(const char *p) {
    if (isam_is_view_index_path(p) || isam_is_store_index_path(p))
        return xstrdup(p);
    return isam_store_index_path_for(p);
}

static uint32_t get_be32(const unsigned char **p) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v = (v << 8) | *(*p)++;
    return v;
}

static uint64_t get_be64(const unsigned char **p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 8) | *(*p)++;
    return v;
}

static int aux_array_item_size(uint8_t type) {
    switch (type) {
    case 'c': case 'C': return 1;
    case 's': case 'S': return 2;
    case 'i': case 'I': case 'f': return 4;
    case 'd': return 8;
    default: return -1;
    }
}

static uint32_t get_le32_at(const unsigned char *p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static size_t aux_payload_len(char type,
                              const unsigned char *payload,
                              const unsigned char *end,
                              const char tag[3]) {
    switch (type) {
    case 'A': case 'c': case 'C':
        if (payload + 1 > end) die("truncated aux tag '%s'", tag);
        return 1;
    case 's': case 'S':
        if (payload + 2 > end) die("truncated aux tag '%s'", tag);
        return 2;
    case 'i': case 'I': case 'f':
        if (payload + 4 > end) die("truncated aux tag '%s'", tag);
        return 4;
    case 'd':
        if (payload + 8 > end) die("truncated aux tag '%s'", tag);
        return 8;
    case 'Z': case 'H': {
        const unsigned char *p = payload;
        while (p < end && *p) ++p;
        if (p >= end) die("unterminated aux tag '%s'", tag);
        return (size_t)(p - payload) + 1;
    }
    case 'B': {
        if (payload + 5 > end) die("truncated aux array tag '%s'", tag);
        int item_size = aux_array_item_size(payload[0]);
        if (item_size < 0) die("unsupported aux array subtype '%c' for '%s'",
                               payload[0], tag);
        uint32_t n = get_le32_at(payload + 1);
        size_t bytes = 1 + 4 + (size_t)n * (size_t)item_size;
        if (payload + bytes > end) die("truncated aux array tag '%s'", tag);
        return bytes;
    }
    default:
        die("unsupported aux type '%c' for '%s'", type, tag);
    }
}

static void tag_stat_add(tag_stat_t **items,
                         size_t *n,
                         size_t *cap,
                         const char tag[3],
                         int64_t payload_bytes,
                         int64_t total_bytes) {
    uint16_t code = isam_tag_code(tag);
    for (size_t i = 0; i < *n; ++i) {
        if ((*items)[i].tag_code != code)
            continue;
        ++(*items)[i].entries;
        (*items)[i].payload_bytes += payload_bytes;
        (*items)[i].total_bytes += total_bytes;
        return;
    }
    if (*n == *cap) {
        size_t next = *cap ? *cap * 2 : 16;
        *items = xrealloc(*items, next * sizeof **items);
        *cap = next;
    }
    (*items)[*n] = (tag_stat_t){
        .tag_code = code,
        .entries = 1,
        .payload_bytes = payload_bytes,
        .total_bytes = total_bytes
    };
    ++*n;
}

static int tag_stat_cmp(const void *a, const void *b) {
    const tag_stat_t *x = a;
    const tag_stat_t *y = b;
    if (x->total_bytes != y->total_bytes)
        return x->total_bytes < y->total_bytes ? 1 : -1;
    if (x->entries != y->entries)
        return x->entries < y->entries ? 1 : -1;
    return x->tag_code < y->tag_code ? -1 : (x->tag_code > y->tag_code);
}

typedef struct tag_description {
    const char tag[3];
    const char *description;
} tag_description_t;

static const char *tag_description(const char tag[3]) {
    static const tag_description_t descriptions[] = {
        { "BC", "sample barcode" },
        { "CB", "cell barcode" },
        { "CR", "raw cell barcode" },
        { "CY", "cell barcode qualities" },
        { "HP", "haplotype" },
        { "ML", "base modification probabilities" },
        { "MM", "base modifications" },
        { "OQ", "original base qualities" },
        { "PS", "phase set" },
        { "RG", "read group" },
        { "RX", "UMI sequence" },
        { "SA", "supplementary alignments" },
        { "ec", "effective CCS coverage" },
        { "fi", "forward IPD" },
        { "fn", "forward pass count" },
        { "fp", "forward pulse width" },
        { "mc", "mapped concordance" },
        { "np", "number of passes" },
        { "qe", "query end in ZMW" },
        { "qs", "query start in ZMW" },
        { "ri", "reverse IPD" },
        { "rm", "repeated-match trimmed" },
        { "rn", "reverse pass count" },
        { "rp", "reverse pulse width" },
        { "rq", "predicted read accuracy" },
        { "sn", "A/C/G/T signal-to-noise" },
        { "we", "last-base frame" },
        { "ws", "first-base frame" },
        { "zm", "ZMW hole number" }
    };
    for (size_t i = 0; i < sizeof descriptions / sizeof descriptions[0]; ++i) {
        if (strcmp(tag, descriptions[i].tag) == 0)
            return descriptions[i].description;
    }
    return "";
}

static void print_tag_stats(const char *title,
                            tag_stat_t *items,
                            size_t n,
                            int top_n) {
    qsort(items, n, sizeof *items, tag_stat_cmp);
    printf("\n%s\n", title);
    printf("tag\tentries\tpayload_bytes\ttotal_bytes\tdescription\n");
    size_t limit = n;
    if (top_n > 0 && (size_t)top_n < limit)
        limit = (size_t)top_n;
    for (size_t i = 0; i < limit; ++i) {
        char tag[3];
        isam_tag_decode(items[i].tag_code, tag);
        printf("%s\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%s\n",
               tag,
               items[i].entries,
               items[i].payload_bytes,
               items[i].total_bytes,
               tag_description(tag));
    }
}

static void scan_inline_aux(const uint8_t *aux,
                            int aux_len,
                            tag_stat_t **items,
                            size_t *n,
                            size_t *cap,
                            int64_t *field_count) {
    if (aux_len == 0)
        return;
    const unsigned char *p = aux;
    const unsigned char *end = aux + aux_len;
    while (p < end) {
        if (p + 3 > end)
            die("truncated inline aux blob");
        char tag[3] = { (char)p[0], (char)p[1], 0 };
        char type = (char)p[2];
        const unsigned char *payload = p + 3;
        size_t payload_len = aux_payload_len(type, payload, end, tag);
        size_t total = 3 + payload_len;
        tag_stat_add(items, n, cap, tag, (int64_t)payload_len, (int64_t)total);
        ++*field_count;
        p += total;
    }
    if (p != end)
        die("corrupt inline aux blob");
}

static void print_dbi_stat(MDB_txn *txn,
                           MDB_dbi dbi,
                           const char *name,
                           size_t page_size) {
    if (!dbi)
        return;
    MDB_stat st;
    int rc = mdb_stat(txn, dbi, &st);
    if (rc)
        die("mdb_stat(%s): %s", name, mdb_strerror(rc));
    uint64_t pages = (uint64_t)st.ms_branch_pages +
                     (uint64_t)st.ms_leaf_pages +
                     (uint64_t)st.ms_overflow_pages;
    printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
           name,
           (uint64_t)st.ms_entries,
           pages * (uint64_t)page_size,
           (uint64_t)st.ms_branch_pages,
           (uint64_t)st.ms_leaf_pages,
           (uint64_t)st.ms_overflow_pages);
}

static void print_storage_explanation(const char *type) {
    printf("file_layout\n");
    if (strcmp(type, "alignment") == 0) {
        printf("role\talignment_view\n");
        printf("stored_here\tread dictionary, imported alignment records, alignment indexes, view-local annotations/subsets, and import-carried read aux in read_aux_by_id\n");
        printf("stored_elsewhere\tcanonical SEQ/QUAL/header payloads stay in the original store BAM/CRAM; canonical read identity and store-side annotations live in the linked .ssam\n");
        printf("reconstruction\taln_by_id provides flags/position/CIGAR/mate/inline aux; read_aux_by_id adds per-read import aux; the linked store supplies canonical SEQ/QUAL and store annotations when viewing\n");
    } else if (strcmp(type, "store") == 0) {
        printf("role\tcanonical_store_index\n");
        printf("stored_here\tcanonical read identity, qname lookup, placement hints, store-side annotations/subsets\n");
        printf("stored_elsewhere\tfull primary BAM/CRAM records, including SEQ/QUAL and canonical aux, stay in the original store file\n");
        printf("reconstruction\tthis .ssam resolves read_id and payload offsets; canonical SEQ/QUAL and canonical aux are then read from the original BAM/CRAM\n");
    } else {
        printf("role\tunknown\n");
    }
    printf("\nmetric_notes\n");
    printf("lmdb_tables.entries\tnumber of records in that LMDB table\n");
    printf("lmdb_tables.page_bytes\tphysical LMDB pages used by each table, including B-tree overhead and free space inside pages\n");
    printf("key_value_bytes\tlogical bytes of LMDB keys plus stored values; excludes LMDB page overhead\n");
    printf("raw_record_bytes\tuncompressed logical alignment-record value bytes before zstd\n");
    printf("zstd_value_bytes\tcompressed value bytes currently stored for zstd-encoded records\n");
    printf("cigar_bytes_in_aln_by_id\tCIGAR bytes embedded in alignment-view records; each CIGAR op is 4 bytes before record compression\n");
    printf("inline_aux_bytes_in_aln_by_id\talignment-local aux fields embedded directly in alignment records, usually tags such as SA/mc/rm\n");
    printf("read_aux_by_id\tbundled per-read aux overrides/import-carried tags stored once per read, not once per alignment\n");
    printf("read_tag\tper-read typed annotation table used by annotate; not a duplicate of canonical BAM/CRAM aux\n\n");
}

static void print_lmdb_storage(lmdb_backend_t *b) {
    MDB_stat env_st;
    int rc = mdb_env_stat(b->env, &env_st);
    if (rc)
        die("mdb_env_stat: %s", mdb_strerror(rc));

    MDB_txn *txn = NULL;
    rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc)
        die("txn_begin: %s", mdb_strerror(rc));

    printf("lmdb_tables\n");
    printf("table\tentries\tpage_bytes\tbranch_pages\tleaf_pages\toverflow_pages\n");
    print_dbi_stat(txn, b->meta, "meta", env_st.ms_psize);
    print_dbi_stat(txn, b->canonical_by_id, "canonical_by_id", env_st.ms_psize);
    print_dbi_stat(txn, b->canonical_by_qname, "canonical_by_qname", env_st.ms_psize);
    print_dbi_stat(txn, b->read_by_id, "read_by_id", env_st.ms_psize);
    print_dbi_stat(txn, b->read_by_qname, "read_by_qname", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_by_id, "aln_by_id", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_import_by_order, "aln_import_by_order", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_by_coord, "aln_by_coord", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_unmapped, "aln_unmapped", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_by_read, "aln_by_read", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_by_bin, "aln_by_bin", env_st.ms_psize);
    print_dbi_stat(txn, b->aln_checkpoint, "aln_checkpoint", env_st.ms_psize);
    print_dbi_stat(txn, b->subset_by_name, "subset_by_name", env_st.ms_psize);
    print_dbi_stat(txn, b->subset_read, "subset_read", env_st.ms_psize);
    print_dbi_stat(txn, b->read_tag, "read_tag", env_st.ms_psize);
    print_dbi_stat(txn, b->read_aux_by_id, "read_aux_by_id", env_st.ms_psize);
    print_dbi_stat(txn, b->annotation_meta, "annotation_meta", env_st.ms_psize);

    mdb_txn_abort(txn);
}

static void decode_read_aux_value(lmdb_backend_t *b,
                                  const MDB_val *val,
                                  uint8_t **aux,
                                  size_t *aux_len,
                                  int *encoding) {
    const unsigned char *p = val->mv_data;
    const unsigned char *end = p + val->mv_size;
    if (p >= end || *p++ != 1 || p + 1 + 4 > end)
        die("corrupt read_aux_by_id value");
    *encoding = *p++;
    uint32_t raw_len = get_be32(&p);
    if (*encoding == 0) {
        if (p + raw_len > end)
            die("corrupt raw read_aux_by_id value");
        *aux = xmalloc(raw_len ? raw_len : 1);
        if (raw_len)
            memcpy(*aux, p, raw_len);
        *aux_len = raw_len;
        return;
    }
    if (*encoding == 1) {
        if (!b->read_aux_zstd_dict || b->read_aux_zstd_dict_len == 0)
            die("read_aux_by_id value requires missing zstd dictionary");
        *aux = xmalloc(raw_len ? raw_len : 1);
        ZSTD_DCtx *dctx = ZSTD_createDCtx();
        if (!dctx)
            die("zstd create decompression context failed");
        size_t dsz = ZSTD_decompress_usingDict(dctx, *aux, raw_len,
                                               p, (size_t)(end - p),
                                               b->read_aux_zstd_dict,
                                               b->read_aux_zstd_dict_len);
        ZSTD_freeDCtx(dctx);
        if (ZSTD_isError(dsz))
            die("zstd decompress read_aux_by_id failed: %s", ZSTD_getErrorName(dsz));
        if (dsz != raw_len)
            die("zstd read_aux_by_id size mismatch");
        *aux_len = raw_len;
        return;
    }
    die("unsupported read_aux_by_id encoding %d", *encoding);
}

static void scan_read_aux(lmdb_backend_t *b,
                          tag_stat_t **items,
                          size_t *n,
                          size_t *cap,
                          int64_t *entries,
                          int64_t *raw_aux_bytes,
                          int64_t *key_value_bytes,
                          int64_t *fields,
                          int64_t *zstd_entries,
                          int64_t *zstd_value_bytes) {
    if (!b->read_aux_by_id)
        return;
    MDB_txn *txn = NULL;
    MDB_cursor *cursor = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc)
        die("txn_begin: %s", mdb_strerror(rc));
    rc = mdb_cursor_open(txn, b->read_aux_by_id, &cursor);
    if (rc)
        die("cursor(read_aux_by_id): %s", mdb_strerror(rc));

    MDB_val key = { 0, NULL };
    MDB_val val = { 0, NULL };
    rc = mdb_cursor_get(cursor, &key, &val, MDB_FIRST);
    while (rc == 0) {
        if (key.mv_size != 8)
            die("corrupt read_aux_by_id key");
        const unsigned char *kp = key.mv_data;
        (void)get_be64(&kp);

        uint8_t *aux = NULL;
        size_t aux_len = 0;
        int encoding = 0;
        decode_read_aux_value(b, &val, &aux, &aux_len, &encoding);
        if (aux_len > INT_MAX)
            die("read_aux_by_id blob too large to scan");
        int64_t field_count_before = *fields;
        scan_inline_aux(aux, (int)aux_len, items, n, cap, fields);
        if (*fields == field_count_before && aux_len != 0)
            die("read_aux_by_id blob contains no fields");
        free(aux);

        ++*entries;
        *raw_aux_bytes += (int64_t)aux_len;
        *key_value_bytes += (int64_t)(key.mv_size + val.mv_size);
        if (encoding == 1) {
            ++*zstd_entries;
            *zstd_value_bytes += (int64_t)val.mv_size;
        }

        rc = mdb_cursor_get(cursor, &key, &val, MDB_NEXT);
    }
    if (rc != MDB_NOTFOUND)
        die("cursor_next(read_aux_by_id): %s", mdb_strerror(rc));

    mdb_cursor_close(cursor);
    mdb_txn_abort(txn);
}

static void scan_read_tag(lmdb_backend_t *b,
                          tag_stat_t **items,
                          size_t *n,
                          size_t *cap,
                          int64_t *entries,
                          int64_t *payload_bytes,
                          int64_t *key_value_bytes) {
    if (!b->read_tag)
        return;
    MDB_txn *txn = NULL;
    MDB_cursor *cursor = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc)
        die("txn_begin: %s", mdb_strerror(rc));
    rc = mdb_cursor_open(txn, b->read_tag, &cursor);
    if (rc)
        die("cursor(read_tag): %s", mdb_strerror(rc));

    MDB_val key = { 0, NULL };
    MDB_val val = { 0, NULL };
    rc = mdb_cursor_get(cursor, &key, &val, MDB_FIRST);
    while (rc == 0) {
        if (key.mv_size < 10)
            die("corrupt read_tag key");
        const unsigned char *k = key.mv_data;
        k += 8;
        char tag[3] = { (char)k[0], (char)k[1], 0 };

        const unsigned char *p = val.mv_data;
        const unsigned char *end = p + val.mv_size;
        if (p >= end || *p++ != 1 || p + 1 + 4 > end)
            die("corrupt read_tag value");
        (void)*p++;
        uint32_t raw_len = get_be32(&p);
        if (p + raw_len + 2 > end)
            die("corrupt read_tag value");
        tag_stat_add(items, n, cap, tag, (int64_t)raw_len,
                     (int64_t)(key.mv_size + val.mv_size));
        ++*entries;
        *payload_bytes += (int64_t)raw_len;
        *key_value_bytes += (int64_t)(key.mv_size + val.mv_size);

        rc = mdb_cursor_get(cursor, &key, &val, MDB_NEXT);
    }
    if (rc != MDB_NOTFOUND)
        die("cursor_next(read_tag): %s", mdb_strerror(rc));

    mdb_cursor_close(cursor);
    mdb_txn_abort(txn);
}

static void scan_view_records(lmdb_backend_t *b,
                              tag_stat_t **inline_tags,
                              size_t *n_inline_tags,
                              size_t *cap_inline_tags,
                              int64_t *records,
                              int64_t *cigar_ops,
                              int64_t *cigar_bytes,
                              int64_t *inline_aux_fields,
                              int64_t *inline_aux_bytes) {
    if (!b->aln_by_id || !b->aln_import_by_order)
        return;
    lmdb_view_iter_t *it = lmdb_view_iter_open(b, LMDB_VIEW_ITER_IMPORT,
                                               -1, 0, 0);
    lmdb_view_record_t rec;
    while (lmdb_view_iter_next(it, &rec)) {
        ++*records;
        *cigar_ops += rec.n_cigar;
        *cigar_bytes += (int64_t)rec.n_cigar * 4;
        *inline_aux_bytes += rec.aux_len;
        scan_inline_aux(rec.aux, rec.aux_len, inline_tags, n_inline_tags,
                        cap_inline_tags, inline_aux_fields);
        lmdb_view_record_clear(&rec);
    }
    lmdb_view_iter_close(it);
}

static void scan_view_record_storage(lmdb_backend_t *b,
                                     int64_t *entries,
                                     int64_t *raw_record_bytes,
                                     int64_t *key_value_bytes,
                                     int64_t *zstd_entries,
                                     int64_t *zstd_value_bytes) {
    if (!b->aln_by_id)
        return;
    MDB_txn *txn = NULL;
    MDB_cursor *cursor = NULL;
    int rc = mdb_txn_begin(b->env, NULL, MDB_RDONLY, &txn);
    if (rc)
        die("txn_begin: %s", mdb_strerror(rc));
    rc = mdb_cursor_open(txn, b->aln_by_id, &cursor);
    if (rc)
        die("cursor(aln_by_id): %s", mdb_strerror(rc));

    MDB_val key = { 0, NULL };
    MDB_val val = { 0, NULL };
    rc = mdb_cursor_get(cursor, &key, &val, MDB_FIRST);
    while (rc == 0) {
        if (key.mv_size != 8)
            die("corrupt aln_by_id key");
        const unsigned char *p = val.mv_data;
        const unsigned char *end = p + val.mv_size;
        if (p >= end)
            die("empty aln_by_id value");
        if (*p == 1) {
            *raw_record_bytes += (int64_t)val.mv_size;
        } else if (*p == 2) {
            if (p + 1 + 1 + 4 > end)
                die("corrupt compressed aln_by_id value");
            ++p;
            int encoding = *p++;
            uint32_t raw_len = get_be32(&p);
            *raw_record_bytes += (int64_t)raw_len;
            if (encoding == 1) {
                ++*zstd_entries;
                *zstd_value_bytes += (int64_t)val.mv_size;
            } else if (encoding != 0) {
                die("unsupported aln_by_id encoding %d", encoding);
            }
        } else {
            die("unsupported aln_by_id value version %u", (unsigned)*p);
        }
        ++*entries;
        *key_value_bytes += (int64_t)(key.mv_size + val.mv_size);
        rc = mdb_cursor_get(cursor, &key, &val, MDB_NEXT);
    }
    if (rc != MDB_NOTFOUND)
        die("cursor_next(aln_by_id): %s", mdb_strerror(rc));

    mdb_cursor_close(cursor);
    mdb_txn_abort(txn);
}

int cmd_store_stats(int argc, char **argv) {
    int top_n = 0;
    int storage_only = 0;
    static struct option longopts[] = {
        {"top",          required_argument, 0, 1000},
        {"storage-only", no_argument,       0, 1001},
        {"help",         no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000:
            top_n = atoi(optarg);
            if (top_n < 0)
                die("--top must be non-negative");
            break;
        case 1001:
            storage_only = 1;
            break;
        case 'h':
            usage();
            return 0;
        default:
            usage();
            return 2;
        }
    }
    if (optind >= argc) {
        usage();
        return 2;
    }

    char *isam_path = resolve_isam_path(argv[optind]);
    lmdb_backend_t *lmdb = lmdb_backend_open(isam_path, 1, 0);
    char *type = lmdb_backend_meta_get(lmdb, "isam_type");
    if (!type)
        die("'%s' has no isam_type meta", isam_path);

    printf("isam_path\t%s\n", isam_path);
    printf("isam_type\t%s\n\n", type);
    print_storage_explanation(type);
    print_lmdb_storage(lmdb);

    if (storage_only) {
        free(type);
        lmdb_backend_close(lmdb);
        free(isam_path);
        return 0;
    }

    tag_stat_t *read_tags = NULL;
    size_t n_read_tags = 0, cap_read_tags = 0;
    int64_t read_tag_entries = 0, read_tag_payload_bytes = 0;
    int64_t read_tag_key_value_bytes = 0;
    scan_read_tag(lmdb, &read_tags, &n_read_tags, &cap_read_tags,
                  &read_tag_entries, &read_tag_payload_bytes,
                  &read_tag_key_value_bytes);

    tag_stat_t *read_aux_tags = NULL;
    size_t n_read_aux_tags = 0, cap_read_aux_tags = 0;
    int64_t read_aux_entries = 0, read_aux_raw_bytes = 0;
    int64_t read_aux_key_value_bytes = 0, read_aux_fields = 0;
    int64_t read_aux_zstd_entries = 0, read_aux_zstd_value_bytes = 0;
    scan_read_aux(lmdb, &read_aux_tags, &n_read_aux_tags, &cap_read_aux_tags,
                  &read_aux_entries, &read_aux_raw_bytes,
                  &read_aux_key_value_bytes, &read_aux_fields,
                  &read_aux_zstd_entries, &read_aux_zstd_value_bytes);

    tag_stat_t *inline_tags = NULL;
    size_t n_inline_tags = 0, cap_inline_tags = 0;
    int64_t records = 0, cigar_ops = 0, cigar_bytes = 0;
    int64_t inline_aux_fields = 0, inline_aux_bytes = 0;
    int64_t view_record_entries = 0, view_record_raw_bytes = 0;
    int64_t view_record_key_value_bytes = 0;
    int64_t view_record_zstd_entries = 0, view_record_zstd_value_bytes = 0;
    if (strcmp(type, "alignment") == 0) {
        scan_view_record_storage(lmdb,
                                 &view_record_entries,
                                 &view_record_raw_bytes,
                                 &view_record_key_value_bytes,
                                 &view_record_zstd_entries,
                                 &view_record_zstd_value_bytes);
        scan_view_records(lmdb, &inline_tags, &n_inline_tags, &cap_inline_tags,
                          &records, &cigar_ops, &cigar_bytes,
                          &inline_aux_fields, &inline_aux_bytes);
        printf("\nalignment_payloads\n");
        printf("records\t%" PRId64 "\n", records);
        printf("aln_by_id_entries\t%" PRId64 "\n", view_record_entries);
        printf("aln_by_id_raw_record_bytes\t%" PRId64 "\n", view_record_raw_bytes);
        printf("aln_by_id_key_value_bytes\t%" PRId64 "\n", view_record_key_value_bytes);
        printf("aln_by_id_zstd_dict_bytes\t%zu\n", lmdb->view_record_zstd_dict_len);
        printf("aln_by_id_zstd_entries\t%" PRId64 "\n", view_record_zstd_entries);
        printf("aln_by_id_zstd_value_bytes\t%" PRId64 "\n", view_record_zstd_value_bytes);
        printf("cigar_ops\t%" PRId64 "\n", cigar_ops);
        printf("cigar_bytes_in_aln_by_id\t%" PRId64 "\n", cigar_bytes);
        printf("inline_aux_fields_in_aln_by_id\t%" PRId64 "\n", inline_aux_fields);
        printf("inline_aux_bytes_in_aln_by_id\t%" PRId64 "\n", inline_aux_bytes);
    }

    printf("\nread_tag_payloads\n");
    printf("entries\t%" PRId64 "\n", read_tag_entries);
    printf("payload_bytes\t%" PRId64 "\n", read_tag_payload_bytes);
    printf("key_value_bytes\t%" PRId64 "\n", read_tag_key_value_bytes);

    printf("\nread_aux_by_id_payloads\n");
    printf("entries\t%" PRId64 "\n", read_aux_entries);
    printf("fields\t%" PRId64 "\n", read_aux_fields);
    printf("raw_aux_bytes\t%" PRId64 "\n", read_aux_raw_bytes);
    printf("key_value_bytes\t%" PRId64 "\n", read_aux_key_value_bytes);
    printf("zstd_dict_bytes\t%zu\n", lmdb->read_aux_zstd_dict_len);
    printf("zstd_entries\t%" PRId64 "\n", read_aux_zstd_entries);
    printf("zstd_value_bytes\t%" PRId64 "\n", read_aux_zstd_value_bytes);

    if (n_read_tags)
        print_tag_stats("read_tag_by_tag", read_tags, n_read_tags, top_n);
    if (n_read_aux_tags)
        print_tag_stats("read_aux_by_id_by_tag", read_aux_tags, n_read_aux_tags, top_n);
    if (n_inline_tags)
        print_tag_stats("inline_aux_by_tag", inline_tags, n_inline_tags, top_n);

    free(inline_tags);
    free(read_aux_tags);
    free(read_tags);
    free(type);
    lmdb_backend_close(lmdb);
    free(isam_path);
    return 0;
}
