# libisam C API

The C API is htslib-like: records are returned as real `bam1_t` values and
headers as `sam_hdr_t` values.

Include:

```c
#include <isam/isam.h>
```

## Opening

```c
isamFile *isam_open(const char *path, const char *mode);
isamFile *isam_open_ext(const char *path, const char *mode,
                        const isamOptions *opts);
int isam_close(isamFile *fp);
```

Read mode (`"r"`) is currently supported.

`isamOptions`:

```c
typedef struct isamOptions {
    const char *store_path_override;
    const char *store_index_path_override;
    const char *reference_fasta;
    int load_store_annotations;
    int load_view_annotations;
    int tag_conflict_policy;
    size_t payload_cache_bytes;
    int n_threads;
} isamOptions;
```

`reference_fasta` is used for CRAM payload decoding. `payload_cache_bytes`
controls the canonical payload cache; zero selects the default. `n_threads`
is passed to htslib for the canonical store handle.

## Sequential Read

```c
sam_hdr_t *hdr = isam_hdr_read(fp);
bam1_t *b = bam_init1();

while (isam_read1(fp, hdr, b) == ISAM_OK) {
    /* inspect b with htslib helpers */
}
```

Returned records are reconstructed from:

- alignment-view fields
- canonical QNAME/SEQ/QUAL payload
- alignment aux tags stored in the view

The view CIGAR is authoritative. Hard clipping is reconciled with canonical
SEQ/QUAL during reconstruction; incompatible lengths report an error.

## Index and Iterators

```c
isamIndex *idx = isam_index_load("sample.hg38.isam");
isamItr *it = isam_itr_querys(idx, hdr, "chr1:1000-2000");
isamItr *by_name = isam_itr_query_name(idx, "read001");
isamItr *subset = isam_itr_query_subset(idx, "hifi");
isamItr *region_subset =
    isam_itr_query_region_subset(idx, hdr, "chr1:1000-2000", "hifi");
```

Fetch records:

```c
while (isam_itr_next(fp, it, b) == ISAM_OK) {
    /* use b */
}
```

Destroy iterators and indexes:

```c
isam_itr_destroy(it);
isam_index_destroy(idx);
```

Region iterators require an interval index. Build one with:

```bash
isamtools sort --interval sample.hg38.isam
```

LMDB views support sequential reads, region iterators, QNAME iterators, named
subset iterators, and region-filtered subset iterators through this API.

## Errors

```c
const char *isam_strerror(int errcode);
int isam_errno(isamFile *fp);
```

Common codes:

- `ISAM_OK`
- `ISAM_EOF`
- `ISAM_ERR_IO`
- `ISAM_ERR_BACKEND`
- `ISAM_ERR_MISSING_QNAME`
- `ISAM_ERR_STALE_INDEX`
- `ISAM_ERR_NO_INDEX`
- `ISAM_ERR_UNSUPPORTED`
