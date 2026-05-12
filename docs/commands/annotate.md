# isamtools annotate

Add read-level tags or named QNAME subsets to a store index or alignment view.

```bash
isamtools annotate [options] <isam-or-store>
```

## Options

- `--from FILE`: input TSV or QNAME list.
- `--from-bam FILE`: in subset mode, stream QNAMEs from SAM/BAM/CRAM records.
- `--subset NAME`: subset mode. Without this, tag mode is used.
- `--description TEXT`: optional subset description.
- `--allow-missing`: skip QNAMEs absent from the target `.isam`.
- `-h, --help`: show command help.

## Tag Format

```text
qname<TAB>TAG:TYPE:VALUE<TAB>TAG:TYPE:VALUE
```

## Subset Format

One QNAME per line.

With `--from-bam`, QNAMEs are read directly from alignment records. Duplicate
QNAMEs are harmless because subset membership is keyed by read id.

## Examples

```bash
isamtools annotate sample.hg38.isam --from tags.tsv
isamtools annotate sample.hg38.isam --subset hifi --from hifi.qnames.txt
isamtools annotate sample.hg38.bam --subset hifi --from-bam sample.hifi.hg38.bam
```
