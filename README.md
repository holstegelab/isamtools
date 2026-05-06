# isamtools

`isamtools` stores one canonical BAM/CRAM read-payload store and many compact
`.isam` alignment views. LMDB is the storage backend. The tool is intended for samples
aligned repeatedly to different references or aligner settings, where
duplicating SEQ/QUAL and invariant read payload in every BAM/CRAM is wasteful.

The canonical store keeps QNAME, SEQ, QUAL, and canonical-location hints. Each
alignment view stores alignment-specific fields such as FLAG, TID, POS, MAPQ,
CIGAR, mate fields, and aux tags. `isamtools view` reconstructs SAM/BAM records
by combining the view with the canonical payload.

## Core Workflow

```bash
isamtools index sample.chm13.bam

isamtools import sample.hg38.bam \
  --store sample.chm13.bam \
  --aln sample.hg38.isam

isamtools sort --all sample.hg38.isam

isamtools view sample.hg38.isam chr1:1000000-2000000
```

The canonical store index and alignment view are file-backed LMDB `.isam`
databases:

```bash
isamtools index -o sample.chm13.bam.isam sample.chm13.bam

isamtools import sample.hg38.bam \
  --store sample.chm13.bam \
  --store-index sample.chm13.bam.isam \
  --aln sample.hg38.isam

isamtools sort --all sample.hg38.isam
```

## Commands

- `index`: build or rebuild a canonical store index.
- `import`: import an alternate alignment into a standalone `.isam` view.
- `sort`: build logical coordinate/name/interval indexes.
- `view`: reconstruct SAM/BAM records from an alignment view.
- `annotate`: add read tags or named QNAME subsets.
- `rehint`: rewrite stale `store_read_id` values after a store-index rebuild.
- `update`: mutate metadata paths and aliases.

Detailed command docs live in [docs/commands](docs/commands).

## Identity Model

The canonical store index uses `canonical_read.read_id` as the integer identity
for a read. Alignment views mirror that ID in `read_dict.store_read_id`, while
preserving QNAME so the mapping can be repaired.

If the store index is rebuilt and read IDs change, run:

```bash
isamtools rehint sample.hg38.isam
```

Every QNAME in the view must still exist in the target store index. Missing
QNAMEs are reported as errors.

## Hard Clipping Policy

The view CIGAR is authoritative. If a view record has hard clipping and its
CIGAR query length differs from the canonical SEQ length, reconstruction trims
SEQ/QUAL only when the leading/trailing hard clips exactly explain the
difference. Otherwise `view` reports an incompatible CIGAR/payload error.

This is important because BAM requires SEQ length to equal the query-consuming
CIGAR length excluding hard clips.

## Performance Notes

View reconstruction fetches canonical payloads by integer `store_read_id` and
uses an in-memory payload cache. Repeated alignments for the same read, such as
secondary and supplementary records, avoid redundant payload decoding.

For region queries, build the interval index:

```bash
isamtools sort --interval sample.hg38.isam
```

For full coordinate/name scans, build all indexes:

```bash
isamtools sort --all sample.hg38.isam
```

## Build

```bash
cmake -S . -B build \
  -DISAM_LMDB_SOURCE_DIR=/path/to/openldap/libraries/liblmdb
cmake --build build
```

Dependencies:

- htslib
- LMDB source (`mdb.c` and `midl.c`)
- samtools for the shell tests

### LMDB Source

LMDB is statically compiled into `libisam.a` from source:

```bash
cmake -S . -B build-lmdb \
  -DISAM_LMDB_SOURCE_DIR=/path/to/openldap/libraries/liblmdb
cmake --build build-lmdb
```

When network access is explicitly available, CMake can fetch the source:

```bash
cmake -S . -B build-lmdb \
  -DISAM_LMDB_FETCH=ON
```

The fetch path is deliberately opt-in so ordinary builds do not depend on
network access.

## Tests

```bash
tests/functional.sh
tests/smoke.sh
tests/lmdb_store.sh
tests/lmdb_alignment.sh
tests/lmdb_compat.sh
N=50000 tests/profile_small.sh
```

The tests use tiny synthetic or subset fixtures and avoid the large files under
`data/`.

## C API

The public C API returns htslib-native `sam_hdr_t` and `bam1_t` objects. See
[docs/api.md](docs/api.md).
