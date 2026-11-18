# Temporary BAMs for unmodified applications

An ordinary `pysam`, samtools, or containerized application cannot open an iSAM
view without the compile-time htslib compatibility shim. A standard indexed BAM
is the portable bridge. Keep the CRAM/iSAM bundle as the retained data and
materialize a BAM only for the lifetime of a downstream task.

If the application already supports reference-backed CRAM and the canonical
alignment is the desired input, it can read `canon.cram` directly instead.

For worker/cache tuning, direct-export commands and example measurements,
see [optimizing ISAM to BAM performance](performance.md). This guide covers
consumer execution, references, staging and temporary-data cleanup.

## Run a consumer

The helper requires Python 3.9+, `samtools`, and a built `isamtools` binary.
It expects a bundle containing `canon.cram`, `canon.ssam`, and the selected
`VIEW.isam`. Staging all inputs additionally requires `canon.cram.crai`.
The verified raw M5 directory below, or a reference FASTA and its `.fai`,
must be the **canonical CRAM's reference**, also when the selected view contains
CHM13 alignments.

This helper reference is only for decoding the retained payload. A downstream
caller must separately receive the reference matching the **selected alignment
view**: hg38 for an hg38 view, CHM13 for a CHM13 view. For a CHM13 analysis,
the helper may therefore use canonical hg38 while the consumer uses CHM13.
Do not substitute the caller's CHM13 FASTA for the canonical decoding reference.

For coordinate export, use a current build supporting raw references,
independent decoders, selected-payload caching and inline indexing. The values
below are examples; adjust them to the available resources.

The current helper requires its `--scratch` directory to be under
`/scratch-node`. This is a helper implementation constraint; direct
`isamtools view` exports can use other working directories as shown in the
[performance guide](performance.md).

```bash
python3 tools/run_with_bam.py \
  --bundle /data/sample_bundle --view phased_hifi_only \
  --reference-path /references/canonical_hg38-m5 \
  --scratch /scratch-node/WORKDIR \
  --isamtools /path/to/isamtools --samtools /path/to/samtools \
  --threads 2 --decode-threads 4 --cache-mb 4096 \
  --cache-selected --decode-lookahead 1024 --write-index \
  --strategy coordinate --stage all --stage-reference \
  --report /results/sample.temporary-bam.json \
  -- python3 /scripts/analyze.py '{bam}' /results/sample_output
```

Budget CPU for both worker pools and the main thread, and leave RAM beyond the
cache setting for decoding, references and the consumer. Working storage must
fit staged inputs, the expanded BAM/BAI and consumer temporary files. Size
these from your data and measured peak memory rather than a fixed machine preset.

The helper emits a fast-compressed BAM in coordinate order, then runs
`samtools quickcheck` and creates its BAI. With a current reader, optional
`--write-index` creates the BAI during conversion and skips the separate
`samtools index` scan. Successful indexing is a hard
coordinate-order check: invalid output is not silently re-sorted. Counts from
the BAM index must equal the selected iSAM record count before the consumer
starts. Counts include secondary, supplementary, and unmapped alignments;
they are not counts of distinct read names.

The consumer receives a normal BAM path in place of `{bam}`. Its Python code
can therefore use an unmodified `pysam.AlignmentFile(path, 'rb')`, including
indexed `fetch()`. The helper preserves the reconstructed SEQ/QUAL and aux
tags; it does not remove methylation, read-group, PacBio, or phasing tags.

Permanent consumer outputs must be written **outside** the temporary directory.
`{scratch}` expands to the helper's unique work directory, if the consumer needs
temporary space. The command runs directly without a shell and must wait for
its own work to finish; do not start a detached background consumer. Its exit
code is propagated. After success or a handled failure, the helper removes only
its own directory, including the BAM, BAI, and any staged inputs. Input bundles
are read-only. The optional JSON report records timings, counts, and cleanup.
Abrupt process termination or node loss can still require scratch cleanup.
The report's `conversion_children_cpu_seconds` is aggregate producer plus sort
CPU, excluding staging, separate indexing and the consumer. With
`--write-index`, inline indexing is included in conversion time/CPU and
`index_seconds` is zero; `index_included_in_conversion` records this distinction.

The helper also forwards `--cache-selected` and `--decode-lookahead N` to a
current reader when independent coordinate decoders are selected. These
options are explicit: `--cache-selected` adds a metadata planning pass to keep
only requested payloads and release final uses; `--decode-lookahead 4096`
enlarges the rolling metadata window without enlarging its pending-decode
limit. See [view options](commands/view.md) for constraints and accounting.
These options remain explicit; the controls below do not change the helper's
cache, worker or `auto` strategy defaults.

See the [example measurements](performance.md#example-measurements) for
validated whole-sample controls and their scope. Use them to understand the
tradeoffs, then measure representative inputs on your own system.

## Recommended raw M5 references for coordinate reconstruction

Prefer `--reference-path` for the direct coordinate/independent-decoder route.
These files contain exactly the same reference bases as FASTA: one uppercase,
newline-free sequence file per canonical contig, named by its header M5
checksum. HTSlib memory-maps them, avoiding repeated FASTA loading, copying,
line removal and uppercasing during container seeks. Both formats can live on
SSD; the gain comes from reference access as well as placement. This changes
neither alignment coordinates nor read sequences, qualities or tags.

FASTA remains supported through `--reference FASTA` / `view -T FASTA`; the
command-line defaults and the helper's `auto` strategy have not changed.
Prepare raw references using the built-in command;
it uses the bundled HTSlib and does not require Python, pysam or samtools:

```bash
isamtools ref-index \
  --cram /data/sample_bundle/canon.cram \
  -T /references/canonical_hg38.fa \
  -o reference-m5
```

For direct conversion, `isamtools view -T reference.fa --prepare-reference DIR`
prepares a fresh directory and switches this conversion to raw references in
one step. Reuse it afterwards with `--reference-path DIR`. Plain FASTA
coordinate decoding prints a preparation hint and continues with FASTA; no
large cache is created implicitly. Native preparation also creates a missing
FASTA `.fai`. See [ref-index](commands/ref-index.md) for the full interface.
The Python `tools/reference_path.py` preparation utility remains available for
workflows that use its JSON summary or preparation API.

The output directory must be new and its parent must exist. Preparation verifies
every canonical header SQ length and MD5 against uppercase FASTA bases in
bounded chunks. It writes raw, newline-free sequence files using the literal
M5 filename, including its case. Existing directories are refused and never
deleted; a failed preparation removes only its newly created directory. Input
CRAM and FASTA files are unchanged. Directory names cannot contain `:` or `%`,
which have special meaning in HTSlib's reference lookup.

In the consumer command, replace `--reference FASTA` with
`--reference-path reference-m5`; they are mutually
exclusive. The helper checks complete canonical SQ file/length coverage before
decoding. Missing entries fail, without falling back to the original FASTA UR or
an online resolver. Once prepared, this route does not require the original
FASTA/FAI. Treat the verified reference directory as immutable: the helper
does not rehash its entire contents on each invocation. It sets `REF_PATH` and
disables `REF_CACHE` for its child processes only. A canonical CRAM consumer
also works through this route; the downstream caller still needs its own
alignment reference as described above.

`--stage-reference` copies the required M5 files into the helper's private SSD
directory. Without that flag, the helper leaves the supplied directory in
place and does not remove it. Account for preparation/copy time and reference
bytes in the workflow's resources. One four-view regional study verified
3,366 contigs into 3.22 GB of files in 40.62 seconds; regional decoder CPU fell
by roughly 37–50%.

A completed 31 October 2026 whole **phased HiFi** comparison used identical
all-SSD inputs, the same rolling decoder build, eight decoder workers, three
BAM-compression workers, level-1 BAM output and a 16-GiB cache:

| Reference | Conversion | BAM indexing | Total |
| --- | ---: | ---: | ---: |
| Raw M5, first run | 8m57s | 1m19s | 10m16s |
| SSD FASTA | 24m31s | 1m30s | 26m01s |
| Raw M5, repeat | 8m29s | 1m31s | 10m01s |

All three cases produced 2,500,335 alignments, including 9,323 unmapped, in a
39.78-GB BAM. Total/unmapped counts and an independent 818-record indexed
SAM/SEQ/QUAL/tag probe passed. Raw references made conversion 2.74–2.89 times
faster and conversion plus indexing 2.54–2.60 times faster in this same-node
control. Staging took 3m44s; verified reference preparation took 40s and wrote
3.22 GB, both excluded from the table. Prepare once and reuse the immutable
reference directory across samples. These are one sample's live-node results;
measure other views and thread counts rather than assuming the same factor.

An older sequential `store-sort` control found no raw-reference gain (23m02s
with SSD FASTA versus 24m49s with raw M5). That route scans the canonical file
once instead of repeatedly seeking containers. The new recommendation is
supported for direct coordinate reconstruction; it does not assign that gain
to every strategy or change `auto` silently.

### Shared lab reference archive

The shared reference bundle for these migrated samples is published at:

```text
agh_full_snellius:/processed/pacbio_isam/references/GRCh38_full_analysis_set_plus_decoy_hla/v1/
```

It contains `GRCh38_full_analysis_set_plus_decoy_hla.raw-m5.v1.tar.gz`,
`reference-m5.manifest.json` and a `COMPLETED.json` publication marker.
The archive contains `reference-m5/` and its manifest; the manifest records
all 3,366 canonical contig names, lengths and M5 checksums and the original
FASTA SHA-256. Reference files are shared across samples rather than copied
into every sample archive. This is the canonical hg38 decoding reference,
including when reconstructing a CHM13 alignment view.

With the lab's checksum-verifying `dcache_cp` installed, download and extract
inside a compute allocation. A tape-only archive can require recall first:

```bash
dcache_cp \
  agh_full_snellius:/processed/pacbio_isam/references/GRCh38_full_analysis_set_plus_decoy_hla/v1/GRCh38_full_analysis_set_plus_decoy_hla.raw-m5.v1.tar.gz \
  /scratch-node/USER.ALLOCATION/
tar -xzf /scratch-node/USER.ALLOCATION/GRCh38_full_analysis_set_plus_decoy_hla.raw-m5.v1.tar.gz \
  -C /scratch-node/USER.ALLOCATION/
```

Use the extracted `reference-m5` directory with `--reference-path`. Since it
is already on SSD, omit `--stage-reference` to reuse it directly. The helper
removes its own temporary BAM/input directory, leaving that supplied reference
directory available for subsequent samples. `dcache_cp` verifies transfer
Adler-32; publication also verifies every raw sequence's length and MD5 after
packaging, and records the archive SHA-256 in the publication marker.
The published archive is 983,109,483 bytes (3,217,346,917 raw sequence bytes);
its SHA-256 is
`dbcbf00da5b3d1187f9d745fb2b7eae2a837a6f47d574b98879d27422c94a44d`.

## Completed scaling and remaining-cost controls

The subsequent matched whole-view comparison completed with all inputs and
references on SSD, level-1 BAM compression and the same 2,500,335 phased HiFi
records:

| Input | Worker configuration | Conversion | Including indexing |
| --- | --- | ---: | ---: |
| ISAM | 1 decoder + 3 output workers | 41m29s | 43m54s |
| ISAM | 2 decoders + 3 output workers | 29m28s | 30m55s |
| ISAM | 4 decoders + 3 output workers | 15m37s | 16m49s |
| ISAM | 8 decoders + 3 output workers | 8m38s | 9m59s |
| Standalone CRAM | 11 HTSlib workers, first/repeat | 2m28s–4m15s | 3m48s–5m58s |

The standalone CRAM was prepared from exactly the selected phased HiFi view;
it excludes the canonical store's extra non-HiFi alignments. Its preparation
time is separate. These are the available-tool comparison (samtools/HTSlib
1.21 versus the ISAM reader's bundled HTSlib 1.23.1), on a shared live node.
The CRAM repeat shows substantial wall-time variation. CRAM normalizes `=`/`X`
CIGAR spelling to `M`; those controls explicitly compare equivalent CIGARs
while retaining all other operation/length, SEQ/QUAL and typed-tag checks.
ISAM output retains the original per-view CIGAR. Consequently the CRAM-derived
BAM is 39.32 GB versus the original-CIGAR ISAM BAM's 39.78 GB.

Separate same-node remaining-cost controls used eight decoders and raw M5
references. The baseline/repeat conversion took 8m24s–8m33s:

| Change from baseline | Conversion | Interpretation |
| --- | ---: | --- |
| Compressed BAM to `/dev/null` | 7m53s | Removing filesystem writes saves about 7%. |
| Uncompressed BAM to `/dev/null` | 7m23s | Removing writes and compression saves about 13% together. |
| 32-GiB instead of 16-GiB cache | 7m46s | About 8% faster, with roughly 40 versus 24 GiB RSS. |
| 4,096 instead of 512 metadata records of look-ahead | 7m50s | About 7–8% faster; no additional decoder workers. |
| Seven instead of three output workers | 7m50s | About 7–8% faster, using a larger CPU reservation. |

These are individual controls, not an additive or combined speedup. Null-sink
cases are diagnostics, not usable BAM products. Every real BAM passed the
record/unmapped counts and the independent indexed payload/tag probe; temporary
data were removed. The default look-ahead and cache settings remain unchanged.
The 32-GiB cache reduced decodes from 46,291 to 42,411, but substantial repeated
decoding remains. This decoder work is a more promising target than SSD
bandwidth; further gains from compression/output tuning are relatively modest.
Indexing still adds roughly 1–1.5 minutes and is reported separately.

## Input placement and scope

- `--stage indexes` (default for iSAM): copy the selected `.isam` and `.ssam`
  to private SSD files, leaving CRAM on shared storage. Copies do not include
  the original LMDB lockfiles.
- `--stage none`: leave inputs in place. Use this for canonical CRAM or for
  iSAM/SSAM files whose resolved targets are under `/scratch-node`. A local
  directory containing symlinks to shared-storage indexes is rejected before
  counting or opening LMDB. CRAM and reference symlinks to shared storage remain allowed.
- `--stage all`: copy the CRAM, CRAI, selected iSAM view, and SSAM to SSD.
- `--stage-reference`: additionally copy the canonical FASTA/FAI, or the
  required M5 files when using `--reference-path`, to SSD.
- `--strategy auto` (default): use `store-sort` for a whole iSAM view and
  `coordinate` for a region or canonical CRAM. Whole-view mode requires a
  current build with `view --order store`; use `--isamtools` to select it.
- `--strategy coordinate`: reconstruct records directly in indexed
  coordinate order, without another sort.
- `--strategy store-sort`: stream the whole canonical input
  once, reconstruct each read's alignments and coordinate-sort the pipe on SSD.
  Requires a build with `view --order store`. Use only for a whole iSAM view,
  not canonical CRAM or `--region`. A full HiFi placement comparison and
  non-HiFi hg38, HiFi CHM13 and non-HiFi CHM13 public-helper cases have passed.
  Both the producer and sorter must
  succeed before indexing/consumer.
- `--region chr20:10000001-20000000`: emit only the requested region. Use this
  only if the downstream analysis genuinely accepts a regional input.

Reserve SSD for the peak simultaneous input copies, expanded BAM, BAI,
downstream temporary files, and headroom. A 71 GB CRAM produced a 103 GB fast
BAM in one whole-canonical measurement; compressed input size is not a safe
estimate of scratch demand. Input-copy time can outweigh the gain for a single
small query. Whole-view comparisons and repeated-query workloads should guide
whether staging pays off; more threads are not automatically faster for random
CRAM container reads.

The `store-sort` strategy runs two processes concurrently. Reserve roughly
`2 * (threads + 1)` CPU cores, plus memory for CRAM decoding and the sort's
512 MiB per worker: the default `--threads 3` whole-view mode needs about eight
cores, not three. Its peak SSD requirement includes sort runs **and** the
merged BAM, in addition to staged inputs and downstream temporary files. It
does not retain an additional unsorted BAM. More allocated threads are not a
substitute for measuring the decoder/sorter balance on the target node.

One completed whole-HiFi study reconstructed 2,500,335 alignments into a
39.76 GB BAM in 19m51s including sort/index, plus 17.59s for private indexes/
reference staging. Copying the 71.05 GB canonical CRAM to SSD took a further
152.13s and did not improve conversion time (20m17s all-SSD, excluding copies).
For that single task, leaving CRAM on GPFS was preferable; staging can still
pay off across repeated tasks or reduce shared-filesystem traffic. This is
one live-node study, not a guaranteed time or a small-region recommendation.
For a region, default `auto` selects `coordinate`; an explicit
`--strategy coordinate` with `--region` also works.

A separate all-SSD non-HiFi hg38 test reconstructed and indexed 4,893,275
alignments into a 60.64 GB BAM in 34m27s, excluding input copies. An unmodified
indexed pysam consumer passed its total count and an independent 605-record
SEQ/QUAL/tag probe. Initial copying of a three-view bundle/reference took 8m32s
on that node and is reused by the later cases. This is not a non-HiFi GPFS/SSD
comparison. HiFi CHM13 additionally passed 2,372,059 alignments in a 37.88 GB
BAM in 22m43s with a matching independent 145-record probe. Non-HiFi CHM13
passed 4,795,811 records in a 59.09 GB BAM in 35m49s, with a matching independent
210-record probe. All three cases are complete; temporary BAM/BAI/sort files
and the shared private input-copy directory were independently confirmed
removed. The original retained blood/brain/child bundle remains untouched.

LMDB explicitly does not support concurrent processes on different hosts
sharing its database/lockfiles on a remote filesystem. In a multi-node GPFS
benchmark, a reader was observed waiting in its LMDB reader-slot mutex while
its BAM output stopped. Private SSD copies avoid that shared reader table
without disabling locking or changing retained database files. The CRAM itself
and its FASTA may remain on GPFS. This is why index staging is the safe default,
not merely a speed tuning option.

For the canonical alignment, select `--view canonical` and provide
`--expected-records N` from a validated migration/index receipt. For a regional
canonical BAM, `N` must be the count for that **region**, not the whole file.
The helper does not scan the entire CRAM again to discover a count. Canonical
regional queries require CRAI; the default index-staging choice becomes `none`
for canonical CRAM. Never substitute the SSAM unique-payload count for an
alignment count.

## Direct coordinate cache and independent decoders

The coordinate route emits an already coordinate-ordered BAM by following the
iSAM index. It does not invoke `samtools sort`, even when canonical hg38 payload
order differs from the selected CHM13 view. Cache misses affect read/decode
work, not output order. Container descriptors now follow the immutable SSAM
container table, with a constant-time linked LRU and a last-container fast
path; the CLI reconstruction borrows the read-only cached payload rather than
making intermediate full-record copies. Public copy-based API behavior is
unchanged.

Use the [resource-sizing guidance](performance.md#choose-a-worker-and-memory-budget)
with verified raw references for the phased-HiFi route above. The helper's
whole-view `auto` policy remains `store-sort`; select `coordinate` explicitly.
The decoder option is rejected for canonical inputs or store-sort.

Independent decoder handles run at most N distinct container tasks from a
rolling metadata window (default 512 alignments, configured with
`--decode-lookahead`). The main thread refills workers before reconstructing
and writing current records, including across window boundaries;
it waits only for the demanded container, not all speculative work. At most
`2 * N` running or completed-but-unconsumed containers are retained outside
the LRU. They cannot evict a borrowed current payload. Even one worker overlaps
decoding with reconstruction/output. Output stays serial in view-index order. Reserve
approximately `decode_threads + threads + 1` cores, cache memory, per-worker
reference/decoder memory, and temporary decoded containers awaiting admission
to the cache. The cache byte budget is not a hard process-RSS bound: a single
oversized container is retained and in-flight decoder results consume extra
memory. Private SSD LMDB copies remain required for multi-node use. No new
whole-data sidecar, unsorted BAM or sort runs are retained by this route.

`ISAM_FETCH_STATS=1` reports actual container decodes and decoded records.
With prefetching, a high demand-fetch hit rate includes preloaded records; it
does not imply that decoding was free. `decode-stats.seconds` sums task wall
time (overlapping for independent workers), not end-to-end conversion time.
The rolling pipeline holds speculative results outside the cache, so cache
hits exclude those first demanded fetches. `pipeline-stats` additionally records
submitted/completed tasks, demand-wait time, and the peak/bound of pending
containers. Completed-task counts include draining any unused tasks at close.
The old and rolling readers' hit percentages are not directly comparable.

The controlled coordinate benchmark compares the frozen pre-change reader,
4/16 GiB caches, BAM-compression worker counts, and subsequently 1/2/4 independent
container decoders. Regional cases validate every output record's SAM/SEQ/QUAL/
tags and successful BAM indexing without sorting; full cases validate alignment
counts plus indexed probes. Copies and scratch cleanup are automatic and
separately reported. Do not extrapolate regional results to complete samples.
The separate `DecodePipelineSnakefile` compares the frozen batch reader with
the rolling pipeline on four 10 Mb regions, then whole HiFi BAMs, using the
same 16 GiB cache and BAM-compression settings. Existing frozen benchmark and
migration executables are not replaced by this experiment.

## Tests

```bash
TOOL=/path/to/isamtools python3 -m unittest discover \
  -s tests -p test_run_with_bam.py -v
```

The small fixtures test complete tag/read-payload preservation, supplementary
and unmapped counts, both strategies, staging, consumer/producer failure,
count-mismatch rejection, regional store-sort refusal, read-only source files,
and temporary cleanup. M5 fixtures additionally cover both strategies,
canonical decoding with the FASTA hidden, missing-reference rejection and
safe preparation failure/existing-directory refusal. No large private dataset
is needed for the test suite.
