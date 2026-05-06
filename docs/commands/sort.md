# isamtools sort

Build logical indexes for a store index or alignment-view `.isam`.

```bash
isamtools sort [options] <file.isam|store.bam|store.cram>
```

## Options

- `-n`: build name/read indexes.
- `--coord`: build coordinate-clustered mapped/unmapped tables.
- `--interval`: build the R*Tree interval index. This also builds coordinate tables.
- `--all`: build all indexes.
- `-h, --help`: show command help.

## Alignment Views

For alignment views, `sort` does not reorder the input file on disk like
`samtools sort`. It builds SQLite tables and indexes used for coordinate scans,
name scans, and interval overlap queries.

For LMDB alignment views, `sort` populates ordered LMDB DBIs for coordinate and
unmapped scans, marks the existing read-name fanout as available, and builds a
fixed-bin interval DBI for region queries.

## Examples

```bash
isamtools sort sample.hg38.isam
isamtools sort -n sample.hg38.isam
isamtools sort --all sample.hg38.isam
```
