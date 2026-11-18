#include "common.h"
#include "reference.h"

#include <htslib/faidx.h>
#include <htslib/hts_log.h>
#include <htslib/kstring.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

int isam_reference_prepare(sam_hdr_t *header, const char *fasta, const char *directory) {
    if (!directory[0] || strchr(directory, ':') || strchr(directory, '%')) {
        hts_log_error("reference directory must be nonempty and contain neither ':' nor '%%'");
        return -1;
    }
    int nref = sam_hdr_nref(header);
    if (nref <= 0) {
        hts_log_error("canonical header has no reference sequences");
        return -1;
    }
    /* Exclusive creation also refuses symlinks and caller-owned directories. */
    if (mkdir(directory, 0777) < 0) {
        hts_log_error("cannot create new reference directory '%s': %s", directory, strerror(errno));
        return -1;
    }
    faidx_t *fai = NULL;
    FILE *output = NULL;
    hts_md5_context *md5 = NULL;
    char **created = xcalloc((size_t)nref, sizeof *created);
    int count = 0, status = -1, temporary_created = 0;
    uint64_t bytes = 0;
    kstring_t m5 = {0, 0, NULL}, path = {0, 0, NULL}, temporary = {0, 0, NULL};
    fai = fai_load3(fasta, NULL, NULL, FAI_CREATE);
    if (!fai) {
        hts_log_error("cannot load/index reference FASTA '%s'", fasta);
        goto done;
    }
    md5 = hts_md5_init();
    if (!md5) {
        hts_log_error("cannot allocate reference checksum context");
        goto done;
    }
    for (int tid = 0; tid < nref; ++tid) {
        const char *name = sam_hdr_tid2name(header, tid);
        hts_pos_t length = sam_hdr_tid2len(header, tid);
        if (!name || length <= 0 ||
            sam_hdr_find_tag_id(header, "SQ", "SN", name, "M5", &m5) < 0 || m5.l != 32) {
            hts_log_error("canonical SQ requires a name, positive length and M5 checksum");
            goto done;
        }
        for (size_t i = 0; i < m5.l; ++i) {
            if (!isxdigit((unsigned char)m5.s[i])) {
                hts_log_error("invalid canonical M5 for '%s'", name);
                goto done;
            }
        }
        if (faidx_seq_len64(fai, name) != length) {
            hts_log_error("FASTA length differs from canonical reference '%s'", name);
            goto done;
        }
        path.l = temporary.l = 0;
        if (ksprintf(&path, "%s/%s", directory, m5.s) < 0 ||
            ksprintf(&temporary, "%s.tmp", path.s) < 0) {
            hts_log_error("cannot allocate raw-reference path");
            goto done;
        }
        output = fopen(temporary.s, "wbx");
        if (!output) {
            hts_log_error("cannot create raw reference '%s': %s", temporary.s, strerror(errno));
            goto done;
        }
        temporary_created = 1;
        hts_md5_reset(md5);
        const hts_pos_t chunk_size = 8 * 1024 * 1024;
        for (hts_pos_t start = 0; start < length;) {
            hts_pos_t wanted = length - start;
            if (wanted > chunk_size) wanted = chunk_size;
            hts_pos_t fetched = 0;
            char *sequence = faidx_fetch_seq64(fai, name, start, start + wanted - 1, &fetched);
            if (!sequence || fetched != wanted) {
                free(sequence);
                hts_log_error("cannot fetch complete reference chunk for '%s'", name);
                goto done;
            }
            for (hts_pos_t i = 0; i < fetched; ++i)
                sequence[i] = (char)toupper((unsigned char)sequence[i]);
            hts_md5_update(md5, sequence, (unsigned long)fetched);
            size_t written = fwrite(sequence, 1, (size_t)fetched, output);
            free(sequence);
            if (written != (size_t)fetched) {
                hts_log_error("cannot write reference '%s'", temporary.s);
                goto done;
            }
            start += fetched;
        }
        unsigned char digest[16];
        char hex[33];
        hts_md5_final(digest, md5);
        hts_md5_hex(hex, digest);
        if (strcasecmp(hex, m5.s)) {
            hts_log_error("FASTA MD5 differs from canonical reference '%s'", name);
            goto done;
        }
        int close_status = fclose(output);
        output = NULL;
        if (close_status != 0) {
            hts_log_error("cannot finish reference '%s'", temporary.s);
            goto done;
        }
        /* A repeated literal M5 can share the file created for a prior SQ. */
        int duplicate = 0;
        for (int i = 0; i < count; ++i)
            if (!strcmp(created[i], path.s)) duplicate = 1;
        if (duplicate) {
            if (unlink(temporary.s) < 0) {
                hts_log_error("cannot remove duplicate temporary reference '%s'", temporary.s);
                goto done;
            }
        } else {
            if (rename(temporary.s, path.s) < 0) {
                hts_log_error("cannot publish reference '%s': %s", path.s, strerror(errno));
                goto done;
            }
            created[count++] = xstrdup(path.s);
            bytes += (uint64_t)length;
        }
        temporary_created = 0;
        temporary.l = 0;
        temporary.s[0] = 0;
    }
    fprintf(stderr, "isamtools ref-index: verified %d contigs, wrote %d files (%" PRIu64 " bytes) to %s\n",
            nref, count, bytes, directory);
    status = 0;
done:
    if (output) fclose(output);
    if (status != 0) {
        if (temporary_created) unlink(temporary.s);
        for (int i = 0; i < count; ++i) unlink(created[i]);
        if (rmdir(directory) < 0)
            hts_log_error("could not remove failed preparation directory '%s': %s", directory, strerror(errno));
    }
    for (int i = 0; i < count; ++i) free(created[i]);
    free(created);
    free(m5.s);
    free(path.s);
    free(temporary.s);
    if (md5) hts_md5_destroy(md5);
    if (fai) fai_destroy(fai);
    return status;
}

static void usage(void) {
    fputs("Usage: isamtools ref-index --cram canon.cram -T reference.fa -o reference.m5\n"
          "\nCreate verified uppercase, newline-free M5 files for fast CRAM payload decoding.\n"
          "The output directory must be new. All canonical SQ lengths and MD5s are checked.\n"
          "A missing FASTA .fai is created. Reuse the result with view --reference-path DIR.\n"
          "\n  --cram FILE      Canonical CRAM supplying reference names, lengths and M5s.\n"
          "  -T, --reference FILE  Reference FASTA.\n"
          "  -o, --output DIR      New raw-reference directory.\n"
          "  -h, --help            Show this help.\n", stderr);
}

int cmd_ref_index(int argc, char **argv) {
    const char *canonical = NULL, *reference = NULL, *directory = NULL;
    static const struct option options[] = {
        {"cram", required_argument, NULL, 1000},
        {"reference", required_argument, NULL, 'T'},
        {"output", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'}, {NULL, 0, NULL, 0}
    };
    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "T:o:h", options, NULL)) != -1) {
        if (c == 1000) canonical = optarg;
        else if (c == 'T') reference = optarg;
        else if (c == 'o') directory = optarg;
        else if (c == 'h') { usage(); return 0; }
        else { usage(); return 2; }
    }
    if (!canonical || !reference || !directory || optind != argc) { usage(); return 2; }
    samFile *input = sam_open(canonical, "r");
    if (!input) die("cannot open canonical CRAM '%s'", canonical);
    const htsFormat *format = hts_get_format(input);
    if (!format || format->format != cram) die("--cram requires a CRAM file");
    sam_hdr_t *header = sam_hdr_read(input);
    if (!header) die("cannot read canonical CRAM header '%s'", canonical);
    int status = isam_reference_prepare(header, reference, directory);
    sam_hdr_destroy(header);
    if (sam_close(input) < 0) status = -1;
    return status == 0 ? 0 : 1;
}
