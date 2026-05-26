# isamtools

When a sample is aligned to multiple references (e.g. hg38 and CHM13) or
re-aligned under different aligner settings, every resulting BAM file contains
the same read sequences and base qualities duplicated in full. For long-read
data especially, this redundancy is substantial. `isamtools` eliminates it.

**One canonical BAM/CRAM. Many lightweight alignment views.**

Existing htslib-based tools (assemblers, variant callers, …) can read `.isam`
files directly via a thin compatibility shim — in most cases with no source
changes at all. See [docs/htslib-compat-migration.md](docs/htslib-compat-migration.md).

`isamtools index` builds a compact LMDB store index (`.ssam`) over the first
(canonical) BAM/CRAM, recording QNAME, placement hints, and a BGZF virtual
offset into the BAM for each read. The original BAM/CRAM is kept and remains
the source of SEQ/QUAL — nothing is re-encoded in LMDB.

Every subsequent alignment is imported with `isamtools import`, which stores
only the alignment-specific fields (FLAG, POS, CIGAR, MAPQ, aux tags, …) in a
second `.isam` view file. That secondary BAM/CRAM can then be deleted.
`isamtools view` reconstructs full SAM/BAM records on the fly by reading
alignment data from the view and fetching the canonical payload from the
original BAM/CRAM via the stored offset.

## Quick start

```bash
# Build once (downloads and statically links htslib + LMDB, and vendors bzip2)
make

# 1. Index the canonical store (keeps sample.chm13.bam; builds sample.chm13.ssam)
isamtools index sample.chm13.bam


# 2. Import a second alignment — only alignment fields are stored, not SEQ/QUAL
isamtools import sample.hg38.bam --store sample.chm13.bam
# creates sample.hg38.isam. sample.hg38.bam can now be deleted; 

# 3. Build coordinate / name / interval indexes on the view
isamtools sort --all sample.hg38.isam

# 4. Query as you would a BAM
isamtools view sample.hg38.isam chr1:1000000-2000000
```

Output is standard SAM/BAM, compatible with any downstream htslib-based tool.

## Annotations and subsets

### Adding per-read tags

Canonical tags that should be shared across imported views can be selected at
index time. Source-specific read-wide tags such as haplotype assignment `HP`
and phase set `PS` are carried by imports and stored once per read.

Share canonical tags at index time:

```bash
isamtools index sample.chm13.bam --share-tags RX,MM,ML
```

Or annotate an existing store or view after the fact with a TSV (`QNAME\tTAG:TYPE:VALUE`):

```bash
# tags.tsv example:
# m64367e_220930_224746/18876851/ccs  HP:i:1  PS:i:100200
# m64367e_220930_224746/97978955/ccs  HP:i:2  PS:i:100200
isamtools annotate sample.chm13.bam --from tags.tsv
```

Use `audit-tags` to preview how tags in your BAM will be classified before
indexing:

```bash
isamtools audit-tags sample.chm13.bam
```

### Named subsets

A subset is a named list of QNAMEs stored inside an alignment `.isam` or store
`.ssam`. `view --subset` checks the alignment view first and then falls back to
the canonical store, so store-level subsets can be reused across multiple
views.

```bash
# Store a list of HiFi read names as the "hifi" subset
isamtools annotate sample.chm13.bam --subset hifi --from hifi.qnames.txt

# Query only HiFi reads from any alignment view
isamtools view --subset hifi sample.hg38.isam chr1:1000000-2000000

# Inspect the canonical store index
isamtools store_stats sample.chm13.ssam
```

Store-level subsets cover all linked views; view-local subsets still work too.

## When to use isamtools

- A sample is aligned to two or more references.
- A sample is re-aligned after parameter tuning, and you want to keep both
  alignments without storing two full copies of the read data.
- You want to attach per-read annotations (tags, QNAME subsets) that are
  shared across all alignments of the same reads.

A single `isamtools import` step replaces a secondary BAM with a view that
stores only the alignment-specific fields — no SEQ/QUAL, no duplicate payload.
The canonical BAM/CRAM is kept as the single source of read data.

## Commands

| Command | Purpose |
|---------|---------|
| `index` | Build or rebuild a canonical store index from a BAM/CRAM. |
| `import` | Import an alternate alignment as a lightweight `.isam` view. |
| `check` | Verify BAM/CRAM/SAM content against a store or view. |
| `sort` | Build coordinate, name, and/or interval indexes on a view. |
| `view` | Reconstruct SAM/BAM records from a view (streaming or region query). |
| `query` | Select reads and emit multi-source TSV fields. |
| `annotate` | Add per-read tags or named QNAME subsets. |
| `store_stats` | Report LMDB table storage and logical payload counts. |
| `audit-tags` | Report how aux tags in a file will be classified (shared, carried, alignment, or unknown). |
| `rehint` | Repair `store_read_id` references after a store index rebuild. |
| `update` | Edit metadata paths and aliases inside a `.ssam` store or `.isam` view. |

**Aux tag classification.** `index --share-tags TAG1,TAG2,...` adds tags to the
canonical BAM/CRAM reconstruction allow-list, while
`index --no-share-tags TAG1,TAG2,...` removes tags from that allow-list.
`import --carry-tags`, `--keep-alignment-tags`, `--drop-tags`,
and `--drop-source-tags` tune how source tags are stored for one imported view.
Use `audit-tags` to preview the classification before indexing or importing.

Detailed per-command documentation: [docs/commands/](docs/commands/).  
Storage layout and internal design: [docs/data-structures.md](docs/data-structures.md).

## Using `.isam` files with existing htslib programs

`libisam` ships a thin compatibility shim that intercepts a small set of
htslib entry points (`hts_open`, `sam_read1`, `sam_itr_querys`, …) and
dispatches to the ISAM backend when the path ends in `.isam`. In the simplest
case (Path A), you add one compile flag and swap the link libraries — **no
source changes required**.

See [docs/htslib-compat-migration.md](docs/htslib-compat-migration.md) for the
three migration paths (compile-time `-include`, single header swap, or explicit
`htsisam_*` calls) and a worked example for the `otter` assembler.

## Build

```bash
make
make packages      # fetch/build vendored third-party static libs only
```

Downloads htslib, LMDB, and bzip2 into `build/vendor`, builds static
libraries, and links `build/isamtools`. The binary has no runtime dependency on
system-installed htslib, LMDB, or libbz2.

Requirements: `cc`, `make`, `curl` or `wget`, `tar`, and development headers
for `zlib` and `liblzma`. `samtools` is needed only to run the shell tests.

```bash
make test          # run functional tests
make clean         # remove build artifacts (keeps downloaded sources)
make distclean     # remove entire build tree including vendor sources
```

Overrides:

```bash
make BUILD_DIR=out
make HTSLIB_VERSION=1.22.2
```

## Tests

```bash
tests/functional.sh
tests/smoke.sh
tests/lmdb_store.sh
tests/lmdb_alignment.sh
tests/lmdb_compat.sh
```

Tests use small fixtures under `test_data/`. Override the binary with
`TOOL=/path/to/isamtools` if needed.

For private, machine-local validation and profiling against the large files in
`data/`, use the opt-in scripts:

```bash
N=10000 tests/private_data_validate.sh
N=100000 tests/profile_private_data.sh
SYNTHETIC=0 CANON_BAM=data/sample.chm13.bam ALN_BAM=data/sample.hg38.bam \
  REPEATS=5 KEEP_WORK=1 tests/profile_comprehensive.sh
```

## C API

The public C API returns htslib-native `sam_hdr_t` and `bam1_t` objects, so it
slots into any existing htslib-based pipeline. See [docs/api.md](docs/api.md).

## Identity and repair

Each read is assigned a stable integer `read_id` in the canonical store.
Alignment views reference that ID as `store_read_id`. If the store is
recreated or re-indexed and IDs change, repair all dependent views with:

```bash
isamtools index sample.hg38.new.cram
isamtools rehint sample.hg38.isam --store sample.hg38.new.cram
```
