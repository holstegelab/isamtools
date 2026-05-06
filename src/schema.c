#include "schema.h"

/*
 * SQL schema for store and alignment-view .isam files.
 *
 * Identity model:
 *   - canonical_read.read_id is the single source of truth for read identity.
 *   - The view's read_dict mirrors that: PK = store_read_id, plus qname
 *     so the mapping can be rebuilt by qname if the store .isam is rebuilt.
 *   - Every annotation table (subset_read, read_tag) keys on the integer id
 *     directly, on both sides.
 */

static const char *STORE_DDL =
    "CREATE TABLE IF NOT EXISTS meta ("
    "  key TEXT PRIMARY KEY,"
    "  value TEXT NOT NULL"
    ");"

    "CREATE TABLE IF NOT EXISTS ref_contig ("
    "  tid INTEGER PRIMARY KEY,"
    "  name TEXT NOT NULL,"
    "  length INTEGER NOT NULL,"
    "  md5 TEXT"
    ");"

    "CREATE TABLE IF NOT EXISTS canonical_read ("
    "  read_id INTEGER PRIMARY KEY,"
    "  qname TEXT NOT NULL UNIQUE,"
    "  qname_hash INTEGER NOT NULL,"
    "  tid INTEGER,"
    "  pos0 INTEGER,"
    "  end0 INTEGER,"
    "  flag INTEGER NOT NULL,"
    "  is_unmapped INTEGER NOT NULL,"
    "  location_kind TEXT NOT NULL,"
    "  bam_voffset INTEGER,"
    "  cram_container_offset INTEGER,"
    "  cram_slice_offset INTEGER,"
    "  cram_container_no INTEGER,"
    "  cram_slice_no INTEGER,"
    "  record_ordinal INTEGER,"
    "  logical_valid INTEGER NOT NULL DEFAULT 1,"
    "  physical_hint_valid INTEGER NOT NULL DEFAULT 1"
    ");"

    "CREATE TABLE IF NOT EXISTS subset ("
    "  subset_id INTEGER PRIMARY KEY,"
    "  name TEXT NOT NULL UNIQUE,"
    "  description TEXT,"
    "  created_at TEXT"
    ");"

    "CREATE TABLE IF NOT EXISTS subset_read ("
    "  subset_id INTEGER NOT NULL REFERENCES subset(subset_id),"
    "  read_id INTEGER NOT NULL REFERENCES canonical_read(read_id),"
    "  PRIMARY KEY (subset_id, read_id)"
    ");"

    "CREATE TABLE IF NOT EXISTS read_tag ("
    "  read_id INTEGER NOT NULL REFERENCES canonical_read(read_id),"
    "  tag TEXT NOT NULL,"
    "  type TEXT NOT NULL,"
    "  value BLOB NOT NULL,"
    "  source TEXT,"
    "  PRIMARY KEY (read_id, tag)"
    ");"

    "CREATE TABLE IF NOT EXISTS annotation_meta ("
    "  annotation_id INTEGER PRIMARY KEY,"
    "  kind TEXT NOT NULL,"
    "  source TEXT,"
    "  command_line TEXT,"
    "  created_at TEXT"
    ");";

static const char *STORE_INDEXES_DDL =
    "CREATE INDEX IF NOT EXISTS canonical_by_qname "
    "  ON canonical_read(qname_hash, qname);"
    "CREATE INDEX IF NOT EXISTS canonical_by_coord "
    "  ON canonical_read(tid, pos0, end0);";

/* Alignment-view DDL.
 * read_dict mirrors store identity: PK is store_read_id, the qname is kept
 * so the mapping can be re-resolved if the store .isam is rebuilt.
 */
static const char *VIEW_DDL =
    "CREATE TABLE IF NOT EXISTS meta ("
    "  key TEXT PRIMARY KEY,"
    "  value TEXT NOT NULL"
    ");"

    "CREATE TABLE IF NOT EXISTS ref_contig ("
    "  tid INTEGER PRIMARY KEY,"
    "  name TEXT NOT NULL,"
    "  length INTEGER NOT NULL,"
    "  md5 TEXT"
    ");"

    "CREATE TABLE IF NOT EXISTS read_dict ("
    "  store_read_id INTEGER PRIMARY KEY,"
    "  qname TEXT NOT NULL UNIQUE,"
    "  qname_hash INTEGER NOT NULL"
    ");"
    "CREATE INDEX IF NOT EXISTS read_dict_by_qname "
    "  ON read_dict(qname_hash, qname);"

    "CREATE TABLE IF NOT EXISTS aln_record_import ("
    "  aln_id INTEGER PRIMARY KEY,"
    "  store_read_id INTEGER NOT NULL,"
    "  flag INTEGER NOT NULL,"
    "  tid INTEGER,"
    "  pos0 INTEGER,"
    "  end0 INTEGER,"
    "  mapq INTEGER,"
    "  n_cigar INTEGER,"
    "  cigar BLOB,"
    "  mate_tid INTEGER,"
    "  mate_pos0 INTEGER,"
    "  tlen INTEGER,"
    "  aux_codec TEXT NOT NULL DEFAULT 'raw',"
    "  aux BLOB,"
    "  record_order INTEGER NOT NULL,"
    "  is_unmapped INTEGER NOT NULL"
    ");"

    "CREATE TABLE IF NOT EXISTS subset ("
    "  subset_id INTEGER PRIMARY KEY,"
    "  name TEXT NOT NULL UNIQUE,"
    "  description TEXT,"
    "  created_at TEXT"
    ");"
    "CREATE TABLE IF NOT EXISTS subset_read ("
    "  subset_id INTEGER NOT NULL REFERENCES subset(subset_id),"
    "  store_read_id INTEGER NOT NULL REFERENCES read_dict(store_read_id),"
    "  PRIMARY KEY (subset_id, store_read_id)"
    ");"

    "CREATE TABLE IF NOT EXISTS read_tag ("
    "  store_read_id INTEGER NOT NULL REFERENCES read_dict(store_read_id),"
    "  tag TEXT NOT NULL,"
    "  type TEXT NOT NULL,"
    "  value BLOB NOT NULL,"
    "  source TEXT,"
    "  PRIMARY KEY (store_read_id, tag)"
    ");"

    "CREATE TABLE IF NOT EXISTS aln_tag ("
    "  aln_id INTEGER NOT NULL,"
    "  tag TEXT NOT NULL,"
    "  type TEXT NOT NULL,"
    "  value BLOB NOT NULL,"
    "  source TEXT,"
    "  PRIMARY KEY (aln_id, tag)"
    ");"

    "CREATE TABLE IF NOT EXISTS annotation_meta ("
    "  annotation_id INTEGER PRIMARY KEY,"
    "  kind TEXT NOT NULL,"
    "  source TEXT,"
    "  command_line TEXT,"
    "  created_at TEXT"
    ");";

static const char *VIEW_CLUSTERED_DDL =
    "CREATE TABLE IF NOT EXISTS aln_record ("
    "  tid INTEGER NOT NULL,"
    "  pos0 INTEGER NOT NULL,"
    "  record_order INTEGER NOT NULL,"
    "  aln_id INTEGER NOT NULL,"
    "  store_read_id INTEGER NOT NULL,"
    "  end0 INTEGER NOT NULL,"
    "  flag INTEGER NOT NULL,"
    "  mapq INTEGER,"
    "  n_cigar INTEGER NOT NULL,"
    "  cigar BLOB,"
    "  mate_tid INTEGER,"
    "  mate_pos0 INTEGER,"
    "  tlen INTEGER,"
    "  aux_codec TEXT NOT NULL DEFAULT 'raw',"
    "  aux BLOB,"
    "  PRIMARY KEY (tid, pos0, record_order, aln_id)"
    ") WITHOUT ROWID;"

    "CREATE TABLE IF NOT EXISTS aln_unmapped ("
    "  record_order INTEGER NOT NULL,"
    "  aln_id INTEGER NOT NULL,"
    "  store_read_id INTEGER NOT NULL,"
    "  flag INTEGER NOT NULL,"
    "  mapq INTEGER,"
    "  n_cigar INTEGER,"
    "  cigar BLOB,"
    "  mate_tid INTEGER,"
    "  mate_pos0 INTEGER,"
    "  tlen INTEGER,"
    "  aux_codec TEXT NOT NULL DEFAULT 'raw',"
    "  aux BLOB,"
    "  PRIMARY KEY (record_order, aln_id)"
    ") WITHOUT ROWID;"

    "CREATE TABLE IF NOT EXISTS aln_id_lookup ("
    "  aln_id INTEGER PRIMARY KEY,"
    "  tid INTEGER NOT NULL,"
    "  pos0 INTEGER NOT NULL,"
    "  record_order INTEGER NOT NULL"
    ");";

void schema_apply_store(sqlite3 *db) {
    db_exec(db, STORE_DDL);
}

void schema_apply_view(sqlite3 *db) {
    db_exec(db, VIEW_DDL);
}

void schema_build_store_indexes(sqlite3 *db) {
    db_exec(db, STORE_INDEXES_DDL);
}

void schema_build_view_clustered(sqlite3 *db) {
    db_exec(db, VIEW_CLUSTERED_DDL);
}

void schema_build_view_indexes(sqlite3 *db,
                               int build_name,
                               int build_coord,
                               int build_interval) {
    if (build_name) {
        db_exec(db,
            "CREATE INDEX IF NOT EXISTS aln_import_by_store_read_id "
            "  ON aln_record_import(store_read_id, record_order);");
        /* read_dict's UNIQUE(qname) and read_dict_by_qname already cover
         * qname lookup. The interesting "by-name" index is the one that
         * fans out from a read to all its alignments. */
        db_exec(db,
            "CREATE INDEX IF NOT EXISTS aln_by_store_read_id "
            "  ON aln_record(store_read_id, record_order);");
        db_exec(db,
            "CREATE INDEX IF NOT EXISTS aln_unmapped_by_store_read_id "
            "  ON aln_unmapped(store_read_id, record_order);");
    }
    if (build_coord) {
        /* clustered table is already sorted by (tid,pos0,...). nothing extra. */
        (void)0;
    }
    if (build_interval) {
        db_exec(db,
            "CREATE VIRTUAL TABLE IF NOT EXISTS aln_interval USING rtree_i32("
            "  aln_id, tid_min, tid_max, pos_min, pos_max"
            ");");
    }
}
