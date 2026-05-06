#include "../include/isam/isam.h"
#include "common.h"
#include "view.h"
#include "store.h"
#include "rebuild.h"

#include <htslib/hts.h>

struct isamFile {
    view_t       *v;
    sqlite3_stmt *seq_stmt;     /* sequential read statement (lazy) */
    bam1_t       *canon;        /* reusable canonical buffer */
    payload_cache_t *cache;
    int           errcode;
};

struct isamIndex {
    char    *path;
    sqlite3 *db;
    sam_hdr_t *hdr_ref;         /* not owned */
};

struct isamItr {
    sqlite3_stmt *stmt;
    int           done;
    int           by_qname;     /* informational */
};

static int step_to_bam(isamFile *fp, sqlite3_stmt *q, bam1_t *b) {
    int rc = sqlite3_step(q);
    if (rc == SQLITE_DONE) return ISAM_EOF;
    if (rc != SQLITE_ROW) { fp->errcode = ISAM_ERR_SQLITE; return ISAM_ERR_SQLITE; }

    int64_t store_read_id = sqlite3_column_int64(q, 0);
    int flag    = sqlite3_column_int(q, 2);
    int tid     = sqlite3_column_type(q,3)==SQLITE_NULL ? -1 : sqlite3_column_int(q,3);
    int64_t pos = sqlite3_column_type(q,4)==SQLITE_NULL ? -1 : sqlite3_column_int64(q,4);
    int mapq    = sqlite3_column_type(q,5)==SQLITE_NULL ?  0 : sqlite3_column_int(q,5);
    int n_cigar = sqlite3_column_type(q,6)==SQLITE_NULL ?  0 : sqlite3_column_int(q,6);
    const uint32_t *cigar = (const uint32_t *)sqlite3_column_blob(q, 7);
    int mtid    = sqlite3_column_type(q,8)==SQLITE_NULL ? -1 : sqlite3_column_int(q,8);
    int64_t mpos= sqlite3_column_type(q,9)==SQLITE_NULL ? -1 : sqlite3_column_int64(q,9);
    int64_t tlen= sqlite3_column_int64(q, 10);
    const uint8_t *aux = (const uint8_t *)sqlite3_column_blob(q, 11);
    int aux_len = sqlite3_column_bytes(q, 11);

    store_record_t srec;
    if (!store_lookup_read_id(fp->v->store, store_read_id, &srec)) {
        fp->errcode = ISAM_ERR_MISSING_QNAME;
        return ISAM_ERR_MISSING_QNAME;
    }
    int frc = payload_cache_fetch(fp->cache, fp->v->store, &srec, fp->canon);
    if (frc != 0) {
        store_record_clear(&srec);
        fp->errcode = frc;
        return frc;
    }

    isam_rebuild_bam(b, fp->canon, flag, tid, pos, mapq,
                     n_cigar, cigar, mtid, mpos, tlen, aux, aux_len);
    store_record_clear(&srec);
    return ISAM_OK;
}

isamFile *isam_open(const char *path, const char *mode) {
    return isam_open_ext(path, mode, NULL);
}

isamFile *isam_open_ext(const char *path, const char *mode,
                        const isamOptions *opts) {
    if (!mode || mode[0] != 'r') {
        /* v1: only read mode is exposed via the C API. */
        return NULL;
    }
    isamFile *fp = xcalloc(1, sizeof *fp);
    fp->v = view_open_read(path);
    view_link_store(fp->v,
                    opts ? opts->store_path_override       : NULL,
                    opts ? opts->store_index_path_override : NULL);
    if (opts && opts->reference_fasta &&
        hts_set_fai_filename(fp->v->store->sam, opts->reference_fasta) != 0) {
        view_close(fp->v);
        free(fp);
        return NULL;
    }
    if (opts && opts->n_threads > 1)
        hts_set_threads(fp->v->store->sam, opts->n_threads);
    fp->canon = bam_init1();
    fp->cache = payload_cache_create(opts ? opts->payload_cache_bytes : 0);
    return fp;
}

int isam_close(isamFile *fp) {
    if (!fp) return 0;
    if (fp->seq_stmt) sqlite3_finalize(fp->seq_stmt);
    if (fp->canon)    bam_destroy1(fp->canon);
    if (fp->cache)    payload_cache_destroy(fp->cache);
    if (fp->v)        view_close(fp->v);
    free(fp);
    return 0;
}

sam_hdr_t *isam_hdr_read(isamFile *fp) {
    /* Return a duplicate so caller's sam_hdr_destroy is safe. */
    sam_hdr_t *dup = sam_hdr_dup(fp->v->hdr);
    return dup;
}

int isam_hdr_write(isamFile *fp, const sam_hdr_t *hdr) {
    (void)fp; (void)hdr;
    return ISAM_ERR_UNSUPPORTED;
}

int isam_read1(isamFile *fp, sam_hdr_t *hdr, bam1_t *b) {
    (void)hdr;
    if (!fp->seq_stmt) {
        int64_t has_coord = 0;
        meta_get_i64(fp->v->db, "indexed_coord", &has_coord);
        fp->seq_stmt = db_prep(fp->v->db, has_coord
            ? "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
              "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
              "FROM ("
              "  SELECT store_read_id, flag, tid, pos0, mapq, n_cigar, cigar,"
              "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id,"
              "         0 AS unmapped_sort "
              "  FROM aln_record "
              "  UNION ALL "
              "  SELECT store_read_id, flag, NULL, NULL, mapq, n_cigar, cigar,"
              "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id,"
              "         1 AS unmapped_sort "
              "  FROM aln_unmapped"
              ") r "
              "JOIN read_dict d ON d.store_read_id=r.store_read_id "
              "ORDER BY r.unmapped_sort, r.tid, r.pos0, r.record_order, r.aln_id;"
            : "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
              "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
              "FROM aln_record_import r "
              "JOIN read_dict d ON d.store_read_id=r.store_read_id "
              "ORDER BY r.is_unmapped, r.tid, r.pos0, r.record_order, r.aln_id;");
    }
    return step_to_bam(fp, fp->seq_stmt, b);
}

isamIndex *isam_index_load(const char *path) {
    isamIndex *idx = xcalloc(1, sizeof *idx);
    idx->path = xstrdup(path);
    idx->db   = db_open(path, SQLITE_OPEN_READONLY);
    return idx;
}

void isam_index_destroy(isamIndex *idx) {
    if (!idx) return;
    if (idx->db) db_close(idx->db);
    free(idx->path);
    free(idx);
}

isamItr *isam_itr_queryi(isamIndex *idx, int tid, hts_pos_t beg, hts_pos_t end) {
    int64_t has_iv = 0;
    meta_get_i64(idx->db, "indexed_interval", &has_iv);
    if (!has_iv) return NULL;

    isamItr *it = xcalloc(1, sizeof *it);
    it->stmt = db_prep(idx->db,
        "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
        "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
        "FROM aln_interval i "
        "JOIN aln_id_lookup k ON k.aln_id=i.aln_id "
        "JOIN aln_record    r ON r.tid=k.tid AND r.pos0=k.pos0 "
        "                    AND r.record_order=k.record_order "
        "                    AND r.aln_id=k.aln_id "
        "JOIN read_dict     d ON d.store_read_id=r.store_read_id "
        "WHERE i.tid_min<=? AND i.tid_max>=? "
        "  AND i.pos_min<? AND i.pos_max>? "
        "ORDER BY r.tid, r.pos0, r.record_order, r.aln_id;");
    sqlite3_bind_int  (it->stmt, 1, tid);
    sqlite3_bind_int  (it->stmt, 2, tid);
    sqlite3_bind_int64(it->stmt, 3, (int64_t)end);
    sqlite3_bind_int64(it->stmt, 4, (int64_t)beg);
    return it;
}

isamItr *isam_itr_querys(isamIndex *idx, sam_hdr_t *hdr, const char *region) {
    int tid = -1;
    hts_pos_t beg = 0, end = HTS_POS_MAX;
    const char *tail = sam_parse_region(hdr, region, &tid, &beg, &end, 0);
    if (!tail || tid < 0) return NULL;
    return isam_itr_queryi(idx, tid, beg, end);
}

isamItr *isam_itr_query_name(isamIndex *idx, const char *qname) {
    isamItr *it = xcalloc(1, sizeof *it);
    it->by_qname = 1;
    int64_t has_coord = 0;
    meta_get_i64(idx->db, "indexed_coord", &has_coord);
    it->stmt = db_prep(idx->db, has_coord
        ? "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
          "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
          "FROM read_dict d "
          "JOIN ("
          "  SELECT store_read_id, flag, tid, pos0, mapq, n_cigar, cigar,"
          "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id "
          "  FROM aln_record "
          "  UNION ALL "
          "  SELECT store_read_id, flag, NULL, NULL, mapq, n_cigar, cigar,"
          "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id "
          "  FROM aln_unmapped"
          ") r ON r.store_read_id=d.store_read_id "
          "WHERE d.qname_hash=? AND d.qname=? "
          "ORDER BY r.record_order, r.aln_id;"
        : "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
          "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
          "FROM read_dict d "
          "JOIN aln_record_import r ON r.store_read_id=d.store_read_id "
          "WHERE d.qname_hash=? AND d.qname=? "
          "ORDER BY r.record_order, r.aln_id;");
    sqlite3_bind_int64(it->stmt, 1, qname_hash(qname));
    sqlite3_bind_text (it->stmt, 2, qname, -1, SQLITE_TRANSIENT);
    return it;
}

isamItr *isam_itr_query_subset(isamIndex *idx, const char *subset_name) {
    if (!subset_name) return NULL;
    isamItr *it = xcalloc(1, sizeof *it);
    int64_t has_coord = 0;
    meta_get_i64(idx->db, "indexed_coord", &has_coord);
    it->stmt = db_prep(idx->db, has_coord
        ? "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
          "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
          "FROM subset s "
          "JOIN subset_read sr ON sr.subset_id=s.subset_id "
          "JOIN read_dict d ON d.store_read_id=sr.store_read_id "
          "JOIN ("
          "  SELECT store_read_id, flag, tid, pos0, mapq, n_cigar, cigar,"
          "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id,"
          "         0 AS unmapped_sort "
          "  FROM aln_record "
          "  UNION ALL "
          "  SELECT store_read_id, flag, NULL, NULL, mapq, n_cigar, cigar,"
          "         mate_tid, mate_pos0, tlen, aux, record_order, aln_id,"
          "         1 AS unmapped_sort "
          "  FROM aln_unmapped"
          ") r ON r.store_read_id=sr.store_read_id "
          "WHERE s.name=? "
          "ORDER BY r.unmapped_sort, r.tid, r.pos0, r.record_order, r.aln_id;"
        : "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
          "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
          "FROM subset s "
          "JOIN subset_read sr ON sr.subset_id=s.subset_id "
          "JOIN read_dict d ON d.store_read_id=sr.store_read_id "
          "JOIN aln_record_import r ON r.store_read_id=sr.store_read_id "
          "WHERE s.name=? "
          "ORDER BY r.is_unmapped, r.tid, r.pos0, r.record_order, r.aln_id;");
    sqlite3_bind_text(it->stmt, 1, subset_name, -1, SQLITE_TRANSIENT);
    return it;
}

isamItr *isam_itr_query_region_subset(isamIndex *idx,
                                      sam_hdr_t *hdr,
                                      const char *region,
                                      const char *subset_name) {
    if (!subset_name) return isam_itr_querys(idx, hdr, region);
    int64_t has_iv = 0;
    meta_get_i64(idx->db, "indexed_interval", &has_iv);
    if (!has_iv) return NULL;

    int tid = -1;
    hts_pos_t beg = 0, end = HTS_POS_MAX;
    const char *tail = sam_parse_region(hdr, region, &tid, &beg, &end, 0);
    if (!tail || tid < 0) return NULL;

    isamItr *it = xcalloc(1, sizeof *it);
    it->stmt = db_prep(idx->db,
        "SELECT r.store_read_id, d.qname, r.flag, r.tid, r.pos0, r.mapq, r.n_cigar, r.cigar,"
        "       r.mate_tid, r.mate_pos0, r.tlen, r.aux "
        "FROM subset s "
        "JOIN subset_read sr ON sr.subset_id=s.subset_id "
        "JOIN read_dict d ON d.store_read_id=sr.store_read_id "
        "JOIN aln_interval i ON 1=1 "
        "JOIN aln_id_lookup k ON k.aln_id=i.aln_id "
        "JOIN aln_record r ON r.tid=k.tid AND r.pos0=k.pos0 "
        "                 AND r.record_order=k.record_order "
        "                 AND r.aln_id=k.aln_id "
        "                 AND r.store_read_id=sr.store_read_id "
        "WHERE s.name=? "
        "  AND i.tid_min<=? AND i.tid_max>=? "
        "  AND i.pos_min<? AND i.pos_max>? "
        "ORDER BY r.tid, r.pos0, r.record_order, r.aln_id;");
    sqlite3_bind_text (it->stmt, 1, subset_name, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int  (it->stmt, 2, tid);
    sqlite3_bind_int  (it->stmt, 3, tid);
    sqlite3_bind_int64(it->stmt, 4, (int64_t)end);
    sqlite3_bind_int64(it->stmt, 5, (int64_t)beg);
    return it;
}

int isam_itr_next(isamFile *fp, isamItr *itr, bam1_t *b) {
    if (!itr || itr->done) return ISAM_EOF;
    int rc = step_to_bam(fp, itr->stmt, b);
    if (rc == ISAM_EOF) itr->done = 1;
    return rc;
}

void isam_itr_destroy(isamItr *itr) {
    if (!itr) return;
    if (itr->stmt) sqlite3_finalize(itr->stmt);
    free(itr);
}

const char *isam_strerror(int errcode) {
    switch (errcode) {
    case ISAM_OK:                  return "ok";
    case ISAM_EOF:                 return "end of stream";
    case ISAM_ERR_IO:              return "i/o error";
    case ISAM_ERR_SQLITE:          return "sqlite error";
    case ISAM_ERR_FORMAT:          return "format error";
    case ISAM_ERR_MISSING_STORE:   return "canonical store not available";
    case ISAM_ERR_MISSING_QNAME:   return "qname missing from store";
    case ISAM_ERR_STALE_INDEX:     return "store/view index is stale";
    case ISAM_ERR_NO_INDEX:        return "required index is missing";
    case ISAM_ERR_UNSUPPORTED:     return "operation not supported";
    default:                       return "unknown error";
    }
}

int isam_errno(isamFile *fp) {
    return fp ? fp->errcode : 0;
}
