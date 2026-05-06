# isamtools annotate

Add read-level tags or named QNAME subsets to a store index or alignment view.

```bash
isamtools annotate [options] <isam-or-store>
```

## Options

- `--from FILE`: input TSV or QNAME list.
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

## Examples

```bash
isamtools annotate sample.hg38.isam --from tags.tsv
isamtools annotate sample.hg38.isam --subset hifi --from hifi.qnames.txt
```
