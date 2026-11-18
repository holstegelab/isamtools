#include "rebuild.h"
#include <time.h>

static const uint8_t complement[16] = {
    0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
};
static const uint8_t aux[] = {'X', 'Y', 'Z', 'o', 'k', 0};

static bam1_t *canonical(int length, int rotation, int flag) {
    const char *alphabet = "=ACMGRSVTWYHKDBN";
    char *sequence = xmalloc((size_t)length + 1);
    char *quality = xmalloc((size_t)length + 1);
    for (int i = 0; i < length; ++i) {
        sequence[i] = alphabet[(i + rotation) & 15];
        quality[i] = i % 64;
    }
    if (length) quality[length - 1] = (char)255;
    uint32_t cigar = bam_cigar_gen(length, BAM_CMATCH);
    bam1_t *record = bam_init1();
    if (bam_set1(record, 5, "read1", flag | (length ? 0 : BAM_FUNMAP),
                 length ? 0 : -1, length ? 10 : -1, 60,
                 length ? 1 : 0, length ? &cigar : NULL,
                 -1, -1, 0, length, sequence, quality, 0) < 0)
        die("cannot build synthetic canonical record");
    /* Padding is not a base, even if an input record contains nonzero bits. */
    if (length & 1) bam_get_seq(record)[length / 2] |= 0xf;
    free(sequence);
    free(quality);
    return record;
}

static int rebuild_and_verify(bam1_t *out, const bam1_t *source,
                              int start, int length, int reverse) {
    int leading = reverse ? source->core.l_qseq - start - length : start;
    int trailing = reverse ? start : source->core.l_qseq - start - length;
    uint32_t cigar[3];
    int n_cigar = 0;
    if (leading) cigar[n_cigar++] = bam_cigar_gen(leading, BAM_CHARD_CLIP);
    if (length) cigar[n_cigar++] = bam_cigar_gen(length, BAM_CMATCH);
    if (trailing) cigar[n_cigar++] = bam_cigar_gen(trailing, BAM_CHARD_CLIP);
    int flag = source->core.flag ^ (reverse ? BAM_FREVERSE : 0);
    isam_rebuild_bam(out, source, flag, 2, 17, 31, n_cigar, cigar,
                     3, 19, 21, aux, sizeof aux);
    if (out->core.l_qseq != length || out->core.flag != flag ||
        out->core.tid != 2 || out->core.pos != 17 || out->core.qual != 31 ||
        out->core.mtid != 3 || out->core.mpos != 19 || out->core.isize != 21 ||
        out->core.n_cigar != n_cigar || strcmp(bam_get_qname(out), "read1") ||
        memcmp(bam_get_cigar(out), cigar, (size_t)n_cigar * 4) ||
        bam_get_l_aux(out) != sizeof aux || memcmp(bam_get_aux(out), aux, sizeof aux))
        return 0;
    for (int i = 0; i < length; ++i) {
        int position = reverse ? start + length - i - 1 : start + i;
        int base = bam_seqi(bam_get_seq(source), position);
        if (reverse) base = complement[base];
        if (bam_seqi(bam_get_seq(out), i) != base ||
            bam_get_qual(out)[i] != bam_get_qual(source)[position])
            return 0;
    }
    return !(length & 1) || !(bam_get_seq(out)[length / 2] & 0xf);
}

static int verify(void) {
    uint64_t cases = 0;
    bam1_t *out = bam_init1();
    for (int n = 0; n <= 33; ++n) {
        for (int rotation = 0; rotation < 16; ++rotation) {
            for (int orientation = 0; orientation < 2; ++orientation) {
                bam1_t *source = canonical(n, rotation, orientation ? BAM_FREVERSE : 0);
                for (int length = n ? 1 : 0; length <= n; ++length) {
                    for (int start = 0; start <= n - length; ++start) {
                        for (int reverse = 0; reverse < 2; ++reverse) {
                            if (!rebuild_and_verify(out, source, start, length, reverse)) {
                                fprintf(stderr, "rebuild mismatch: n=%d rotation=%d orientation=%d start=%d length=%d reverse=%d\n",
                                        n, rotation, orientation, start, length, reverse);
                                bam_destroy1(source);
                                bam_destroy1(out);
                                return 1;
                            }
                            ++cases;
                        }
                    }
                }
                bam_destroy1(source);
            }
        }
    }
    bam_destroy1(out);
    printf("{\"verified_cases\":%llu}\n", (unsigned long long)cases);
    return 0;
}

static int benchmark(int iterations, int length) {
    if (iterations < 1 || length < 1) die("positive benchmark iterations/length required");
    bam1_t *source = canonical(length, 0, 0);
    bam1_t *out = bam_init1();
    uint32_t cigar = bam_cigar_gen(length, BAM_CMATCH);
    struct timespec begin, end;
    uint64_t checksum = 0;
    clock_gettime(CLOCK_MONOTONIC, &begin);
    for (int i = 0; i < iterations; ++i) {
        isam_rebuild_bam(out, source, 0, 2, 17, 31, 1, &cigar,
                         3, 19, 21, aux, sizeof aux);
        checksum += bam_seqi(bam_get_seq(out), i % length);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    double seconds = end.tv_sec - begin.tv_sec + (end.tv_nsec - begin.tv_nsec) / 1e9;
    printf("{\"iterations\":%d,\"length\":%d,\"seconds\":%.9f,\"checksum\":%llu}\n",
           iterations, length, seconds, (unsigned long long)checksum);
    bam_destroy1(out);
    bam_destroy1(source);
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 1) return verify();
    if (argc == 4 && strcmp(argv[1], "--benchmark") == 0)
        return benchmark(atoi(argv[2]), atoi(argv[3]));
    die("usage: rebuild_sequences [--benchmark ITERATIONS LENGTH]");
}
