#ifndef ISAM_REFERENCE_H
#define ISAM_REFERENCE_H

#include <htslib/sam.h>

/* Create a new raw-M5 directory, verifying every canonical SQ length/MD5.
 * On failure remove only files created by this call; existing paths are refused. */
int isam_reference_prepare(sam_hdr_t *header, const char *fasta, const char *directory);

#endif
