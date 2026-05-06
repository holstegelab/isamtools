# isamtools import

Import an alternate alignment into a compact alignment-view `.isam`.

```bash
isamtools import [options] <input.sam|bam|cram|->
```

Required options:

- `--store PATH`: canonical BAM/CRAM store.
- `--aln PATH`: output alignment-view `.isam`.

Optional:

- `--store-index PATH`: store index path. Defaults to `<store>.isam`.
- `--name NAME`: logical alias for the view.
- `--allow-missing`: skip input records absent from the canonical store.
- `--write-index`: run `sort --all` after import.
- `-h, --help`: show command help.

## Behavior

Each input record is resolved by QNAME against the store index. By default,
missing QNAMEs are errors. Imported records preserve view-specific alignment
fields and aux tags; SEQ/QUAL come from the canonical store during `view`.

`import` writes a file-backed LMDB `.isam` containing typed read-dictionary and
alignment-record DBIs. LMDB creates a sibling lock file named `<path>-lock`
while the environment is in use.

## Examples

```bash
isamtools import sample.hg38.bam \
  --store sample.chm13.bam \
  --aln sample.hg38.isam

minimap2 -ax map-ont hg38.fa reads.fq.gz |
  isamtools import - --store sample.chm13.bam --aln sample.hg38.isam
```
