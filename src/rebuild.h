#ifndef ISAM_REBUILD_H
#define ISAM_REBUILD_H

#include "common.h"

void isam_rebuild_bam(bam1_t *out, const bam1_t *canonical,
                      int flag, int tid, int64_t pos, int mapq,
                      int n_cigar, const uint32_t *cigar,
                      int mtid, int64_t mpos, int64_t isize,
                      const uint8_t *aux, int aux_len);

#endif
