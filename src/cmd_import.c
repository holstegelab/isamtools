#include "common.h"
#include "store.h"
#include "view.h"
#include "../include/isam/isam.h"

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <getopt.h>

static void usage(void) {
    fputs(
        "Usage: isamtools import [options] <input.sam|bam|cram|-> \\\n"
        "         --store <store.bam|cram> --aln <out.isam>\n"
        "\n"
        "Options:\n"
        "  --store PATH         Canonical store (BAM/CRAM).\n"
        "  --store-index PATH   Store .isam (default <store>.isam).\n"
        "  --aln PATH           Output alignment-view .isam.\n"
        "  --name NAME          Logical alias for the alignment view.\n"
        "  --allow-missing      Skip records whose QNAME isn't in the store.\n"
        "  --write-index        Run `sort --all` after import.\n"
        "  -h, --help           Show this help.\n",
        stderr);
}

extern int cmd_sort_main(int argc, char **argv);  /* for --write-index */

int cmd_import(int argc, char **argv) {
    const char *store_path = NULL;
    const char *store_idx  = NULL;
    const char *aln_path   = NULL;
    const char *name_alias = NULL;
    int allow_missing = 0;
    int write_index   = 0;

    static struct option longopts[] = {
        {"store",         required_argument, 0, 1000},
        {"store-index",   required_argument, 0, 1001},
        {"aln",           required_argument, 0, 1002},
        {"name",          required_argument, 0, 1003},
        {"allow-missing", no_argument,       0, 1004},
        {"write-index",   no_argument,       0, 1005},
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
        case 1004: allow_missing = 1;   break;
        case 1005: write_index   = 1;   break;
        case 'h':  usage(); return 0;
        default:   usage(); return 2;
        }
    }
    if (optind >= argc) { usage(); return 2; }
    if (!store_path) die("--store is required");
    if (!aln_path)   die("--aln is required");
    const char *input_path = argv[optind];

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
    sam_hdr_t *in_hdr = sam_hdr_read(in);
    if (!in_hdr) die("cannot read header from input '%s'", input_path);

    /* Open / create output view. */
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
    char *ts = now_iso8601(); view_meta_set(v, "created_at", ts); free(ts);
    view_meta_set(v, "indexed_coord",    "0");
    view_meta_set(v, "indexed_name",     "0");
    view_meta_set(v, "indexed_interval", "0");
    view_meta_set(v, "has_import_table", "1");

    view_set_header(v, in_hdr);

    lmdb_backend_begin(v->lmdb);
    bam1_t *b = bam_init1();
    int64_t n_in = 0, n_kept = 0, n_missing = 0;
    int64_t max_ref_span = 0;
    int64_t max_query_len = 0;
    int64_t batch = 0;
    const int64_t BATCH_SIZE = 50000;

    while (1) {
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
                die("QNAME '%s' not found in canonical store '%s' "
                    "(use --allow-missing to skip)", qname, store_path);
            ++n_missing;
            ++n_in;
            continue;
        }

        int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
        int n_cigar = b->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(b);
        const uint8_t  *aux   = bam_get_aux(b);
        int aux_len = bam_get_l_aux(b);
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
        store_record_clear(&srec);

        ++n_kept;
        ++n_in;
        ++batch;
        if (batch >= BATCH_SIZE) {
            lmdb_backend_commit(v->lmdb);
            lmdb_backend_begin(v->lmdb);
            batch = 0;
            fprintf(stderr, "isamtools import: %lld in / %lld kept\r",
                    (long long)n_in, (long long)n_kept);
        }
    }
    lmdb_backend_commit(v->lmdb);
    bam_destroy1(b);
    sam_hdr_destroy(in_hdr);
    sam_close(in);

    fprintf(stderr, "isamtools import: %lld in, %lld kept, %lld missing\n",
            (long long)n_in, (long long)n_kept, (long long)n_missing);

    view_meta_set_i64(v, "view_max_ref_span", max_ref_span);
    view_meta_set_i64(v, "view_max_query_len", max_query_len);

    view_close(v);
    store_close(store);
    free(store_idx_alloc);

    if (write_index) {
        char *args[] = { "sort", "--all", (char *)aln_path, NULL };
        return cmd_sort_main(3, args);
    }
    return 0;
}
