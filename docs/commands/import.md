# isamtools import

Import an alternate alignment into a compact alignment-view `.isam`.

```bash
isamtools import [options] <input.sam|bam|cram|->
```

Required options:

- `--store PATH`: canonical BAM/CRAM store.

Optional:

- `--aln PATH`: output alignment-view `.isam`. Defaults to `<input-root>.isam`.
- `--store-index PATH`: store index path. Defaults to `<store-root>.ssam`.
- `--name NAME`: logical alias for the view.
- `--carry-tags TAGS`: store comma-separated source tags once per read.
- `--keep-alignment-tags TAGS`: keep comma-separated source tags inline per alignment.
- `--drop-tags TAGS`: suppress tags from both source and canonical/store reconstruction.
- `--drop-source-tags TAGS`: ignore matching imported source tags only.
- `--check`: after writing the view, verify input records match reconstructed output.
- `--strict-tag-policy`: fail on stored-once tag conflicts.
- `--verbose-tag-policy`: print the detailed tag policy summary.
- `--allow-missing`: skip input records absent from the canonical store.
- `--write-index`: run `sort --all` after import.
- `-@, --threads N`: set a total worker budget shared between input
  BAM/CRAM decompression and ordered background compression of view records.
  With `N > 1`, one worker compresses records and the remaining workers decode
  the input; `N = 1` keeps both operations synchronous.
- `--fast-verify`: compare every shared tag with the compact canonical
  fingerprint in a current `.ssam`, avoiding canonical BAM/CRAM payload reads.
  With an older `.ssam`, primary records use adaptive payload sampling while
  secondary and supplementary records are compared exactly.
- `-h, --help`: show command help.

## Behavior

Each input record is resolved by QNAME against the store index. By default,
missing QNAMEs are errors. Imported records preserve view-specific alignment
fields and aux tags; SEQ/QUAL come from the canonical store during `view`.

The canonical store index keeps one payload record per QNAME. It prefers a
primary alignment and otherwise falls back to a non-hard-clipped supplementary
alignment. QNAMEs represented only by secondary or hard-clipped supplementary
records are reported as missing unless `--allow-missing` is used.

The default store-index path replaces a trailing `.bam`, `.cram`, or `.sam`
with `.ssam`. For example, `sample.chm13.bam` uses `sample.chm13.ssam`.

Relative `--store` and `--store-index` input arguments use the caller's working
directory and are persisted as absolute paths. This prevents an output view in
another directory from rebasing those input paths. To create a movable bundle,
use `update` afterward to set paths relative to the view's own directory.

If a previous `import` failed and left an incomplete output `.isam`, rerunning
the same command will remove that incomplete file and recreate it.

`import` writes a file-backed LMDB `.isam` containing typed read-dictionary and
alignment-record DBIs. LMDB creates a sibling lock file named `<path>-lock`
while the environment is in use.

Tag policy precedence is deterministic: `--drop-tags` suppresses final output,
then values from the imported alignment file override canonical-backed values,
then missing canonical-backed tags can be filled from the canonical BAM/CRAM
payload.
`--drop-source-tags` removes only the imported copy, so canonical shared
fallback can still happen.

`HP` and `PS` are alignment-local by default. A read can align in more than one
phase set, and haplotype labels are only meaningful within their associated
phase set. Use `--carry-tags HP,PS` only when the source is known to define one
assignment for the complete read. If an explicitly carried tag has conflicting
values, the first observed value is retained and the conflict is reported;
with `--strict-tag-policy`, that conflict is an error.

Alignment-local tags such as `SA`, `mc`, and `rm` are stored on each alignment
record. The import summary separates these from tags stored once per read in
`read_aux_by_id`.

`import` prints the intended tag-policy behavior before scanning, then prints an
observed tag summary at the end. The observed summary only lists tags that were
actually seen.

`--check` is a validation flag for normal imports. It still writes the `.isam`
view, then streams the input and verifies that every input record is
represented by the reconstructed view. Input aux tags must survive with the
same values; output may contain additional policy-derived tags such as
canonical shared fallback.

`--fast-verify` is intended for trusted read-level tag policies such as PacBio
`MM`, `ML`, and `RG`. An index built by the current `isamtools index` stores a
64-bit fingerprint for each canonical shared tag. Import compares shared fields
on primary, secondary and supplementary records with those fingerprints and
stores a canonical-absent value as a read-level override, without decoding the
canonical payload. A genuine within-read disagreement is reported as a
stored-once conflict. Because fingerprints are probabilistic identifiers, omit
`--fast-verify` when exact byte comparison is required; the default path still
fetches and compares the canonical values directly. Older store indexes remain
readable and automatically use the previous adaptive-sampling path.

## Examples

```bash
isamtools import sample.hg38.bam \
  --store sample.chm13.bam

minimap2 -ax map-ont hg38.fa reads.fq.gz |
  isamtools import - --store sample.chm13.bam --aln sample.hg38.isam

isamtools import sample.hg38.bam \
  --store sample.chm13.bam \
  --threads 4 \
  --fast-verify \
  --drop-source-tags RX

isamtools import --check sample.hg38.bam \
  --store sample.chm13.bam \
  --aln sample.hg38.isam
```
