# isamtools view

Reconstruct SAM/BAM from an alignment-view `.isam`.

```bash
isamtools view [options] <view.isam> [region]
```

For whole-view performance, use the [performance guide](../performance.md):
verified raw references, private SSD indexes, explicit worker/cache budgets,
selected-payload caching and inline BAI. It explains how to choose settings
for your hardware and measure their effect.

## Options

- `-b, --bam`: output BAM.
- `-C, --cram`: output CRAM; requires `-T`.
- `-1, --fast`: fast BAM compression and default output to BAM.
- `-u, --uncompressed`: uncompressed BAM output.
- `-O, --output-fmt sam|bam|cram[,level=N|fast|uncompressed]`: choose output format.
- `--output-fmt-option level=N|fast|uncompressed`: output compression option.
- `-h, --header, --with-header`: include the header in SAM output.
- `-H, --header-only`: print only the header.
- `--no-header`: suppress SAM header output. This is the default for SAM, matching `samtools view`.
- `--store PATH`: override canonical store path.
- `--store-index PATH`: override store index path.
- `-T, --reference PATH`: reference FASTA for CRAM payload decoding and CRAM output.
- `--reference-path DIR`: alternative input-only reference lookup using raw,
  uppercase, newline-free files named by the canonical header's M5 checksums.
  These files must have been checksum-verified when populated. Every canonical
  SQ entry must have a matching file of the correct length; missing entries
  fail rather than falling back to the original UR or an online resolver.
  Uses HTSlib's mmap-backed local `REF_PATH` route and disables `REF_CACHE` for
  this process. Cannot be combined with `-T`; CRAM output still requires `-T`.
  Filenames preserve literal header M5 case; directory names cannot contain
  `:` or `%`, which have special lookup semantics. Prepare verified files with
  [`isamtools ref-index`](ref-index.md). The [temporary-BAM helper](../temporary-bam.md)
  accepts the same option. Prefer this reference path for direct coordinate
  reconstruction: a complete phased-HiFi control with eight decoders reduced
  conversion from 24m31s (SSD FASTA) to 8m29s–8m57s (SSD raw M5). See the guide
  for measured scope, setup costs and the shared reference archive. FASTA
  remains available through `-T`; reference contents are unchanged.
- `--prepare-reference DIR`: with `-T FASTA`, prepare a **new** raw-M5 directory
  against the canonical CRAM header, then use it for this conversion. Requires
  SAM/BAM record output; count/header-only requests, CRAM output, canonical BAM
  and `--reference-path` cannot be combined with this option. Preparation
  verifies every reference length and MD5. The directory is kept for reuse with
  `--reference-path DIR`. Plain FASTA coordinate decoding instead prints a
  short hint and continues without creating a cache implicitly.
- `-@, --threads N`: htslib output I/O threads and BAM/legacy input I/O threads.
  By default CRAM container-locator payload fetches decode synchronously to avoid
  repeatedly discarding threaded read-ahead on seeks; output compression remains
  threaded. A sequential `--order store` scan also uses threaded CRAM decoding.
- `-L, --targets-file PATH`, `--regions-file PATH`: BED file of target intervals.
- `--decode-threads N`: opt-in independent CRAM container workers, `1..64`.
  Requires an indexed canonical CRAM and coordinate output, either a whole
  coordinate-indexed view or one region (including `*`). A metadata rolling
  look-ahead (default 512, configured with `--decode-lookahead`) supplies
  container tasks; each worker owns its input handle
  and reference state. Workers are refilled while current records are rebuilt
  and written. Only a demanded container is waited for, not an entire decode
  batch. Pending running/completed containers are bounded at `2 * N` and held
  outside the LRU until demanded. Even `N=1` uses a background decoder.
  Only payload reads are parallel: BAM/SAM records are written
  in the original view-index order, without an external sort. `-@` independently
  controls output compression. Not supported with BED, subsets, fetch-pairs,
  count, header-only or `--unmap`; invalid combinations fail explicitly.
- `-M, --use-index`: accepted for `samtools view` compatibility; indexed paths are used automatically.
- `-P, --fetch-pairs`: for region or BED selection, expand selected read names to all alignments for those reads.
- `--order coordinate|name|store`: full-scan order. `store` streams the
  canonical file once and reconstructs each read's view alignments from its
  exact indexed payload locator. Output is explicitly `SO:unsorted` and must
  be coordinate-sorted before indexing. It requires current BAM/CRAM payload
  locators and cannot be combined with regions, BED, subsets, or fetch-pairs.
  Filters and tag reconstruction are unchanged. Coordinate remains the CLI
  default; see [performance tuning](../performance.md) for measured strategy choices.
- `--subset NAME`: emit only reads in a named subset, searching the alignment view first and then the canonical store.
- `-f INT`: require all FLAG bits.
- `-F INT`: exclude records with any FLAG bits.
- `--rf INT`: require at least one FLAG bit.
- `-G INT`: exclude records with all FLAG bits.
- `-q INT`: minimum MAPQ.
- `-m INT`: minimum query bases consumed by the CIGAR.
- `-r STR`: require read group `RG:Z:STR`.
- `-R PATH`: require read group listed in `PATH`.
- `-N PATH`: require QNAME listed in `PATH`.
- `-l STR`: require library `STR`, resolved through the output header `@RG ID -> LB`.
- `-d TAG[:VALUE]`: require an aux tag, optionally with a matching value.
- `-D TAG:PATH`: require an aux tag whose value appears in `PATH`.
- `--subsample F`: keep fraction `F` of reads, from `0` to `1`.
- `--subsample-seed INT`: seed used by `--subsample`.
- `-s INT.FRAC`: samtools-compatible shorthand for `--subsample 0.FRAC --subsample-seed INT`.
- `-x TAGS`, `--remove-tag TAGS`: strip comma-separated aux tags from output.
  `-x ^TAGS` preserves only those tags, matching samtools shorthand.
- `--keep-tag TAGS`: preserve only comma-separated aux tags in output.
- `-p`, `--unmap`: emit records not selected by filters/regions with the
  `UNMAP` flag set instead of dropping them.
- `-c, --count`: count records instead of emitting them.
- `-o, --output PATH`: output path. Defaults to stdout.
- `--decode-lookahead N`: coordinate metadata window, `1..65536` (default
  `512`), requiring `--decode-threads`. This bounds metadata read-ahead rather
  than the number of decoded containers; pending decodes remain bounded at
  `2 * decode_threads`. Larger windows may keep workers busier across changes
  in canonical container order, without changing output order.
- `--cache-selected`: opt-in exact payload request planning, requiring
  `--decode-threads`. An additional metadata-only pass counts each selected
  canonical locator's uses. Decoding retains only requested full payloads, and
  the main thread frees each payload after its final use. This avoids keeping
  unused non-HiFi reads and already-consumed records in the decoded cache.
  Shared reads remain available for all secondary/supplementary alignments;
  SEQ/QUAL and typed tags are reconstructed as before. The plan uses four bytes
  per canonical record in referenced containers plus a small descriptor table.
  No persistent plan or source change is made. `ISAM_FETCH_STATS=1` reports
  planning time/bytes, requests, remaining requests and final-use releases.
  Planning time is part of conversion time. This option is not enabled by
  default; see the [performance guide](../performance.md) for the completed
  phased-HiFi controls, cache costs and worker tuning.
- `--write-index`: build `OUTPUT.bai` while writing a named coordinate BAM,
  using HTSlib's writer index. It avoids reading the completed BAM a second
  time. Whole views or one region are supported; stdout/nonregular outputs,
  SAM/CRAM, store/name order, BED, subsets, fetch-pairs, count/header-only and
  `--unmap` combinations are refused. HTSlib checks coordinate order and BAI
  bounds; index-save and output-close failures return an error. The BAI is
  ready on successful return, so a separate `samtools index` is unnecessary.
- `-S`, `--no-PG`, and `-?` are accepted for `samtools view` compatibility.

LMDB alignment views support full scans, name-ordered scans, filters, counts,
BED filtering, one or more positional regions, and region queries after
`isamtools sort --interval` or `--all`. Named subsets are searched in the
alignment view first and then in the linked canonical store.

An unfiltered whole-view `--count` reads the LMDB table cardinality directly.
It does not scan alignment records or open the canonical BAM/CRAM and store
index. Counts with regions, subsets, or filters still scan the required index
and record fields. Whole-view filtered counts scan the primary alignment table
sequentially rather than dereferencing the coordinate index record by record.

## Store paths

Relative `store_path` and `store_index_path` metadata are resolved from the
directory containing the `.isam`, so an archive bundle can be moved as a unit.
Absolute metadata paths are used unchanged. Explicit `--store` and
`--store-index` arguments retain ordinary command-line semantics: relative
arguments are resolved from the caller's working directory.

For example, metadata `canon.cram` and `canon.ssam` in
`/path/to/bundle/hifi.hg38.isam` refer to those files in
the same bundle directory. Readers built before this path fix need explicit
overrides, or must run from the bundle directory.

## Hard Clipping

Primary canonical reads in the store must not be hard-clipped. Hard clipping is
allowed in alignment views. The view CIGAR is authoritative; if hard clips make
the view query length shorter than the canonical payload, `view` trims
canonical SEQ/QUAL only when the hard clips exactly account for the difference.
Other length mismatches are errors.

## Examples

```bash
isamtools view sample.hg38.isam chr1:1000000-2000000
isamtools view -h sample.hg38.isam chr1:1000000-2000000
isamtools view -b -1 -o out.bam sample.hg38.isam
isamtools view -C -T ref.fa -o out.cram sample.hg38.isam
isamtools view -c -L targets.bed sample.hg38.isam
isamtools view -c -N qnames.txt -r rg1 sample.hg38.isam chr1:1-100000 chr2:1-100000
isamtools view -P sample.hg38.isam chr1:1000000-2000000
isamtools view -c -F 2304 -q 20 sample.hg38.isam
isamtools view --subsample 0.1 --subsample-seed 42 sample.hg38.isam
isamtools view -x MM,ML sample.hg38.isam
isamtools view -p sample.hg38.isam chr1:1000000-2000000
isamtools view --subset hifi sample.hg38.isam chr1:1-100000
```

For a whole coordinate-indexed view, a current CRAM-backed reader can write
the BAM and BAI together. These example values should be tuned to the available
CPU, memory and storage:

```bash
ISAM_SLICE_CACHE_MB=4096 isamtools view --reference-path ref.m5 \
  --store local_canon.cram --store-index local_canon.ssam \
  --decode-threads 4 --decode-lookahead 1024 --cache-selected -@ 2 \
  --write-index -1 -o output.bam local_view.isam
```

For automatic staging, validation and cleanup around an application, use the
[temporary-BAM helper](../temporary-bam.md). Full recipes and resource sizing
are in the [performance guide](../performance.md).

An explicit sequential alternative scans the canonical store and sorts:

```bash
set -o pipefail
isamtools view --order store -u -@ 3 -T ref.fa local_view.isam |
  samtools sort -@ 3 -m 512M -l 1 -T ./sort_tmp -o output.bam -
samtools index output.bam
```

Use private node-local `.isam`/`.ssam` copies for multi-node workloads; see
[temporary BAMs](../temporary-bam.md). This scans the whole canonical input,
so it is not a replacement for efficient small regional queries.
