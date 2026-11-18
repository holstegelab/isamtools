# Optimizing ISAM to BAM performance

Tune reference access, decoding, caching and output separately. The best
settings depend on the selected alignments, available CPU and memory, and
storage throughput. Measure complete exports, including any indexing and input
copies needed by the application.

For CRAM-backed coordinate export, the main options are:

| Area | Option | Purpose |
| --- | --- | --- |
| Reference access | `--reference-path DIR` | Use verified raw genome references to reduce repeated FASTA processing. |
| Decoding | `--decode-threads N` | Decode independent CRAM containers concurrently. |
| Read-ahead | `--decode-lookahead N` | Look further ahead in coordinate metadata to keep decoders supplied. |
| Payload retention | `--cache-selected` | Retain requested payloads and release them after their final use. |
| Cache budget | `ISAM_SLICE_CACHE_MB` | Bound retained decoded payloads in MiB. |
| Output compression | `-@ N`, `-1` | Configure output workers and fast BAM compression. |
| Indexing | `--write-index` | Build the BAI during writing, avoiding a separate BAM scan. |

These options require a current build; check `isamtools view --help` for
availability. The view needs a coordinate index, and the canonical `.ssam`
needs CRAM payload locators. See [view indexing](commands/sort.md) and
[canonical store indexing](commands/index.md) when preparing a new bundle.

## Choose a worker and memory budget

Start with a modest decoder count and increase it until elapsed time stops
improving. Tune output compression workers independently. Both pools share
the same CPU budget, alongside the main reconstruction/writer thread:

```text
configured CPU workers = decode_threads + output_threads + 1
```

Use that total when sizing a job or sharing a machine with other work. It
describes configured concurrency; average CPU use can be lower. Doubling the
number of decoders does not imply twice the throughput once reconstruction,
compression or writing limits progress.

If compression uses most CPU, compare more output workers or a faster
compression setting. Faster compression can produce a larger BAM, so measure
output size and storage traffic as well as elapsed time.

Choose a cache budget that leaves memory for references, decoder buffers,
pending containers, metadata, compression and the downstream application.
The cache setting is not a process memory limit. Measure peak resident memory
on representative input before increasing concurrency or cache size.

Keep the metadata window and decoder count separate: `--decode-lookahead`
defaults to 512 records, while pending decoded containers are bounded at
`2 * decode_threads`. A larger window can expose enough independent work to
keep workers busy, at the cost of more buffered metadata.

## Prepare the canonical genome reference once

```bash
isamtools ref-index --cram sample/canon.cram -T reference.fa -o reference.m5
```

This verifies reference names, lengths and M5 checksums against the canonical
CRAM, then writes uppercase, newline-free sequence files. The output directory
must be new. Reuse the immutable directory for stores using the same canonical
genome reference. See [ref-index](commands/ref-index.md) for details and
one-step `--prepare-reference` use.

Raw references allow HTSlib to memory-map bases, reducing repeated FASTA loading,
copying, line removal and uppercasing during container seeks. FASTA remains
available through `-T`; compare the two on your access pattern if reference
handling is a substantial part of decoding time.

Use the **canonical CRAM's genome reference** for payload decoding. The
application consuming the exported BAM needs the reference matching that
alignment view. These can differ: a CHM13 alignment view may obtain its read
payloads from a canonical store encoded against hg38.

## Choose input placement

Fast local storage can help when decoding revisits many containers. Include
the cost of copying inputs when deciding whether staging pays off, especially
for a single small query. Reuse staged references and inputs across repeated
tasks when practical.

Keep `.isam` and `.ssam` files private to each host in multi-host workloads.
Sharing LMDB lockfiles across hosts on a network filesystem can leave readers
waiting on a process-shared mutex. The canonical CRAM and reference files can
remain shared if their storage performance is sufficient.

Working storage must fit staged inputs, the expanded BAM/BAI, downstream
temporary files and headroom. Estimate these from your data; compressed input
size alone does not predict the required space.

## Export directly in coordinate order

The following values illustrate configuration; adjust them to your CPU and
memory budget and compare with the defaults:

```bash
decode_threads=4
output_threads=2
cache_mb=4096
lookahead=1024

ISAM_SLICE_CACHE_MB="$cache_mb" isamtools view \
  --reference-path reference.m5 \
  --store sample/canon.cram --store-index sample/canon.ssam \
  --order coordinate --decode-threads "$decode_threads" \
  --decode-lookahead "$lookahead" --cache-selected \
  -@ "$output_threads" -1 --write-index \
  -o output.bam sample/alignments.isam

samtools quickcheck -v output.bam
samtools idxstats output.bam
```

`--cache-selected` adds a metadata pass that counts canonical payload uses.
Decoding retains complete requested payloads, and each is released after its
final use. Secondary and supplementary alignments retain access to their
shared payloads. SEQ/QUAL and per-view tag reconstruction are preserved.
Include planning overhead when judging whether this policy helps your query.

Output follows the view's coordinate index even when canonical payload order
differs. `--write-index` creates `output.bam.bai` during writing, so another
`samtools index` pass is unnecessary. It requires a named coordinate BAM.

Independent decoders support a whole view or one region, including unmapped
records. BED, subsets, fetch-pairs and some other modes require a different
execution path; see the [view option constraints](commands/view.md).

For an application that expects a normal indexed BAM,
[the temporary-BAM helper](temporary-bam.md) automates input staging, validation,
consumer execution and cleanup. It accepts the same tuning options, with
`--threads` for output workers and `--cache-mb` for the payload cache. Select
`--strategy coordinate` explicitly: its whole-view `auto` strategy uses
`store-sort`. The helper documentation also describes its current scratch-path
requirements.

## Example measurements

One phased long-read sample produced a 39.78-GB BAM containing 2,500,335
alignments, including 9,323 unmapped. The following controls used raw references,
local SSD input/output, a 16-GiB cache and level-1 BAM compression. Planning and
inline indexing are included; staging and reference preparation are excluded.
Baseline runs bracketed the optimized cases.

| Configuration | Decoders / output workers | Look-ahead | Conversion | Separate indexing | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline, first/repeat | 8 / 3 | 512 | 7m58s–8m49s | 1m26s–1m43s | 9m24s–10m32s |
| Selected cache | 8 / 3 | 512 | 6m19s | 1m27s | 7m45s |
| Selected cache + inline BAI | 8 / 7 | 4096 | 5m19s | Included | 5m19s |
| Selected cache + inline BAI | 16 / 7 | 4096 | 4m31s | Included | 4m31s |

The fastest case used about 55% less elapsed time than the mean baseline.
Container decodes fell from 46,291 to 23,690 with selected caching. Increasing
from eight to 16 decoders saved about 15% elapsed time, illustrating diminishing
returns rather than a universal thread-count recommendation. Both optimized
inline-index cases peaked at about 23 GiB resident memory.

A separate coordinate-decoding control found raw references 2.74–2.89 times
faster than FASTA with both on SSD. Another control combined a larger cache,
more look-ahead and more output workers without selected caching; it roughly
doubled resident memory for a small total-time improvement. These results
motivate measuring each change against its resource cost.

All accepted BAMs passed quickcheck, indexed total/unmapped counts and an
independent 818-record SAM/SEQ/QUAL/typed-tag probe, including phasing,
methylation, read-group and PacBio tags. Regional tests also compared every
record in an interval larger than the metadata ring. Exact runtimes and gains
will vary with input, hardware, competing work and storage.

## Profile your workload

Compare complete exports on the same system with consistent input, output
compression, selection and validation. Change one setting at a time, then test
useful combinations. Repeat a baseline to expose load and filesystem-cache
variation. Small regional timings do not predict whole-sample throughput.

- Record elapsed time, aggregate CPU time and peak RSS. For example, GNU
  `time -v` reports these around an export command.
- Measure staging and reference preparation separately, then include them when
  comparing end-to-end application runs.
- Count indexing consistently. With inline BAI, the helper reports
  `index_included_in_conversion: true` and `index_seconds: 0`.
- Set `ISAM_FETCH_STATS=1` to inspect container decodes, decoded records,
  planning costs, pending bounds and demand waits. A completed selected-payload
  plan should have zero remaining requests.
- `decode-stats.seconds` sums overlapping task times; it is not elapsed export
  time. Prefetching also changes what a cache hit means, so inspect actual
  decode counts alongside hit percentages.
- If demand waits remain high, test more look-ahead or decoders. If waits are
  low, profile reconstruction, compression and writing before adding workers.
  Remeasure the bottleneck after each substantial improvement.
