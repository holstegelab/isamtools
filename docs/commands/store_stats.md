# isamtools store_stats

```bash
isamtools store_stats [options] <file.isam|file.ssam|store.bam|store.cram>
```

Report where an `.isam` alignment view or `.ssam` canonical store spends storage.

## Options

- `--top N`: show only the top N tags per tag table. By default, show all tags.
- `--storage-only`: skip logical scans and report only LMDB page usage.
- `-h, --help`: show command help.

## Output

`store_stats` prints:

- A `file_layout` section explaining which payloads live in this LMDB and
  which are read from the linked canonical BAM/CRAM or store `.ssam`.
- A `metric_notes` section explaining the main byte counters.
- LMDB table page usage and entry counts.
- For alignment views, CIGAR operation count and CIGAR bytes stored in
  `aln_by_id`.
- Raw and zstd-compressed `aln_by_id` view-record storage.
- Inline alignment aux bytes stored in `aln_by_id`.
- Per-read annotation payloads stored in `read_tag`, grouped by tag.
- Bundled import read-level aux payloads stored in `read_aux_by_id`, grouped by
  tag.

The `payload_bytes` columns count only logical aux value payloads. The
`total_bytes` columns include the LMDB key and encoded value bytes, but not LMDB
page slack or B-tree overhead. Use `lmdb_tables.page_bytes` for the physical
page footprint.

## Examples

```bash
isamtools store_stats sample.chm13.ssam
isamtools store_stats --storage-only sample.hg38.isam
isamtools store_stats --top 20 sample.hg38.isam
```
