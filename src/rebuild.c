#include "rebuild.h"

static int cigar_query_len(int n_cigar, const uint32_t *cigar,
                           int *leading_hclip, int *trailing_hclip) {
    int qlen = 0;
    *leading_hclip = 0;
    *trailing_hclip = 0;

    if (n_cigar <= 0 || !cigar) return 0;

    if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)
        *leading_hclip = bam_cigar_oplen(cigar[0]);
    if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CHARD_CLIP)
        *trailing_hclip = bam_cigar_oplen(cigar[n_cigar - 1]);

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (bam_cigar_type(op) & 1) qlen += len;
    }
    return qlen;
}

static int cigar_has_hard_clip(int n_cigar, const uint32_t *cigar) {
    for (int i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)
            return 1;
    }
    return 0;
}

static hts_pos_t cigar_ref_span(int n_cigar, const uint32_t *cigar) {
    hts_pos_t span = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (bam_cigar_type(op) & 2) span += len;
    }
    return span;
}

static uint8_t nt16_complement(uint8_t base) {
    static const uint8_t comp[16] = {
        0, 8, 4, 12, 2, 10, 6, 14,
        1, 9, 5, 13, 3, 11, 7, 15
    };
    return comp[base & 0xf];
}

static void copy_seq_slice(uint8_t *dst, const uint8_t *src,
                           int start, int len, int reverse) {
    memset(dst, 0, (size_t)((len + 1) / 2));
    for (int i = 0; i < len; ++i) {
        int src_i = reverse ? start + len - 1 - i : start + i;
        uint8_t base = bam_seqi(src, src_i);
        bam_set_seqi(dst, i, reverse ? nt16_complement(base) : base);
    }
}

static void copy_qual_slice(uint8_t *dst, const uint8_t *src,
                            int start, int len, int reverse) {
    if (!reverse) {
        memcpy(dst, src + start, (size_t)len);
        return;
    }

    for (int i = 0; i < len; ++i)
        dst[i] = src[start + len - 1 - i];
}

void isam_rebuild_bam(bam1_t *out, const bam1_t *canonical,
                      int flag, int tid, int64_t pos, int mapq,
                      int n_cigar, const uint32_t *cigar,
                      int mtid, int64_t mpos, int64_t isize,
                      const uint8_t *aux, int aux_len) {
    if (cigar_has_hard_clip(canonical->core.n_cigar, bam_get_cigar(canonical)))
        die("canonical store record '%s' has hard clipping; "
            "re-index from an unclipped primary alignment",
            bam_get_qname(canonical));

    int leading_hclip = 0, trailing_hclip = 0;
    int cigar_qlen = cigar_query_len(n_cigar, cigar,
                                     &leading_hclip, &trailing_hclip);
    int canonical_qlen = canonical->core.l_qseq;
    int out_qlen = canonical_qlen;
    int seq_start = 0;
    int reverse = ((flag ^ canonical->core.flag) & BAM_FREVERSE) != 0;

    if (n_cigar > 0 && cigar_qlen != canonical_qlen) {
        if (cigar_qlen >= 0 &&
            leading_hclip + cigar_qlen + trailing_hclip == canonical_qlen) {
            seq_start = reverse ? trailing_hclip : leading_hclip;
            out_qlen = cigar_qlen;
        } else {
            die("view CIGAR query length (%d) is incompatible with canonical "
                "SEQ length (%d)", cigar_qlen, canonical_qlen);
        }
    }

    int l_qname = canonical->core.l_qname;
    size_t seq_bytes = (size_t)((out_qlen + 1) / 2);
    size_t total = (size_t)l_qname + (size_t)n_cigar * 4 + seq_bytes
                 + (size_t)out_qlen + (size_t)aux_len;

    if ((size_t)out->m_data < total) {
        out->data = (uint8_t *)xrealloc(out->data, total);
        out->m_data = total;
    }

    uint8_t *p = out->data;
    memcpy(p, canonical->data, l_qname);
    p += l_qname;

    if (n_cigar > 0) {
        memcpy(p, cigar, (size_t)n_cigar * 4);
        p += (size_t)n_cigar * 4;
    }

    copy_seq_slice(p, bam_get_seq(canonical), seq_start, out_qlen, reverse);
    p += seq_bytes;
    copy_qual_slice(p, bam_get_qual(canonical), seq_start, out_qlen, reverse);
    p += out_qlen;
    if (aux_len > 0) memcpy(p, aux, (size_t)aux_len);

    out->core = canonical->core;
    out->core.flag = flag;
    out->core.tid = tid;
    out->core.pos = pos;
    out->core.qual = mapq;
    out->core.n_cigar = n_cigar;
    out->core.mtid = mtid;
    out->core.mpos = mpos;
    out->core.isize = isize;
    out->core.l_qseq = out_qlen;

    if (tid >= 0 && pos >= 0 && n_cigar > 0) {
        hts_pos_t span = cigar_ref_span(n_cigar, cigar);
        out->core.bin = hts_reg2bin(pos, pos + (span > 0 ? span : 1), 14, 5);
    } else {
        out->core.bin = 0;
    }

    out->l_data = (int)total;
    out->id = 0;
}
