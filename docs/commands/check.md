# isamtools check

Verify BAM/CRAM/SAM content against an `.isam` alignment view or a canonical
store.

```bash
isamtools check [options] <input.sam|bam|cram>
```

## Options

- `--view PATH`: alignment-view `.isam` to compare against.
- `--store PATH`: canonical store BAM/CRAM for store-only checks, or store override for view checks.
- `--store-index PATH`: store `.ssam` path.
- `--batch-records N`: for view checks, collect this many records and group
  canonical payload fetches by CRAM container before comparison (default
  65,536).
- `--allow-missing`: in store-only mode, skip primary records absent from the store.
- `-h, --help`: show command help.

## Behavior

If `--view` is supplied, `check` compares the input records with the
reconstructed view in import order. If `--view` is omitted and `<input-root>.isam`
exists as an alignment view, that view is used automatically.

If no alignment view is available, `--store` is required. In store-only mode,
`check` compares primary input records against the canonical store and reports
secondary/supplementary records as skipped.

Input aux tags must be present with the same values. Reconstructed view output
may contain additional policy-derived tags, for example canonical shared-tag
fallback.

View checks retain input comparison order, but canonical CRAM payloads are
fetched in bounded batches ordered by container offset. This avoids repeatedly
decoding a container when records from several containers are interleaved,
without loading or reordering the complete input. Larger batches improve
container reuse at the cost of holding more decoded input and canonical records
in memory.

## Examples

```bash
isamtools check sample.hg38.bam --view sample.hg38.isam

isamtools check sample.nonhifi.hg38.bam \
  --store sample.hifi_nonhifi.hg38.bam
```
