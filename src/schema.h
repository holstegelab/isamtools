#ifndef ISAM_SCHEMA_H
#define ISAM_SCHEMA_H

#include "common.h"

void schema_apply_store(sqlite3 *db);
void schema_apply_view(sqlite3 *db);

void schema_build_store_indexes(sqlite3 *db);

void schema_build_view_clustered(sqlite3 *db);    /* aln_record (clustered) + aln_unmapped */
void schema_build_view_indexes(sqlite3 *db,
                               int build_name,
                               int build_coord,
                               int build_interval);

#endif
