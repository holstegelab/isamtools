#include "common.h"
#include "backend_lmdb.h"

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>
#include <inttypes.h>

static void usage(void) {
    fputs(
        "Usage: isamtools index [options] <store.bam|store.cram>\n"
        "\n"
        "Build the canonical store index <store>.isam.\n"
        "By default builds the QNAME index (required for `import`).\n"
        "\n"
        "Options:\n"
        "  --name NAME      Logical alias for this store (e.g. chm13).\n"
        "  -o, --output P   Override default output path (<store>.isam).\n"
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

static int index_lmdb(const char *store_path, const char *out_path,
                      const char *name_alias, samFile *sam,
                      sam_hdr_t *hdr, int is_bam,
                      const struct stat *stbuf) {
    struct stat out_st;
    if (stat(out_path, &out_st) == 0)
        die("LMDB output path '%s' already exists; choose a new path", out_path);
    if (errno != ENOENT)
        die("stat '%s': %s", out_path, strerror(errno));

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

    bam1_t *b = bam_init1();
    int64_t n_records = 0;
    int64_t max_read_len = 0;
    int have_voffset = is_bam;

    while (1) {
        int64_t voffset = -1;
        if (have_voffset)
            voffset = bgzf_tell(sam->fp.bgzf);

        int rc = sam_read1(sam, hdr, b);
        if (rc == -1) break;
        if (rc < -1)
            die("error reading '%s' at record %lld",
                store_path, (long long)n_records);
        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
            continue;
        if ((int64_t)b->core.l_qseq > max_read_len)
            max_read_len = b->core.l_qseq;

        int is_unmapped = (b->core.flag & BAM_FUNMAP) ? 1 : 0;
        lmdb_store_record_t rec = {
            .read_id = n_records + 1,
            .qname = bam_get_qname(b),
            .qname_hash = qname_hash(bam_get_qname(b)),
            .tid = b->core.tid,
            .pos0 = (is_unmapped || b->core.pos < 0) ? -1 : (int64_t)b->core.pos,
            .end0 = is_unmapped ? -1 : (int64_t)bam_endpos(b),
            .flag = b->core.flag,
            .is_unmapped = is_unmapped,
            .bam_voffset = (have_voffset && voffset >= 0) ? voffset : -1,
            .has_voffset = (have_voffset && voffset >= 0) ? 1 : 0
        };
        lmdb_store_put_record(lmdb, &rec);

        ++n_records;
        if (n_records % 100000 == 0) {
            lmdb_backend_commit(lmdb);
            lmdb_backend_begin(lmdb);
            fprintf(stderr, "isamtools index: %lld records...\r",
                    (long long)n_records);
        }
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(sam);

    lmdb_store_put_meta_i64(lmdb, "indexed_name", 1);
    lmdb_store_put_meta_i64(lmdb, "indexed_coord", 1);
    lmdb_store_put_meta_i64(lmdb, "store_max_read_len", max_read_len);
    lmdb_backend_commit(lmdb);
    lmdb_backend_close(lmdb);

    fprintf(stderr, "isamtools index: %lld records indexed\n",
            (long long)n_records);
    fprintf(stderr, "isamtools index: wrote LMDB %s\n", out_path);
    return 0;
}

int cmd_index(int argc, char **argv) {
    const char *name_alias = NULL;
    const char *out_path   = NULL;

    static struct option longopts[] = {
        {"name",   required_argument, 0, 1000},
        {"output", required_argument, 0, 'o'},
        {"help",   no_argument,       0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:h", longopts, NULL)) != -1) {
        switch (c) {
        case 1000: name_alias = optarg; break;
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

    samFile *sam = sam_open(store_path, "r");
    if (!sam) die("cannot open '%s'", store_path);
    const htsFormat *fmt = hts_get_format(sam);
    int is_bam = (fmt && fmt->format == bam);
    int is_cram = (fmt && fmt->format == cram);
    if (!is_bam && !is_cram)
        die("unsupported canonical store format (need BAM or CRAM)");

    sam_hdr_t *hdr = sam_hdr_read(sam);
    if (!hdr) die("cannot read header of '%s'", store_path);

    /* metadata */
    struct stat stbuf;
    if (stat(store_path, &stbuf) != 0)
        die("stat '%s': %s", store_path, strerror(errno));

    int rc = index_lmdb(store_path, out_path, name_alias, sam, hdr,
                        is_bam, &stbuf);
    free(out_alloc);
    return rc;
}
