# isamtools update

Mutate metadata in a store index or alignment-view `.isam`.

```bash
isamtools update [options] <isam-or-store>
```

## Options

- `--name NAME`: set logical alias.
- `--store-path PATH`: set canonical store path on an alignment view.
- `--store-index-path PATH`: set store-index path on an alignment view.
- `--reference-name NAME`: set store reference name metadata.
- `--reference-fasta PATH`: set store reference FASTA metadata.
- `--set KEY=VALUE`: set arbitrary metadata.
- `--unset KEY`: remove metadata.
- `-h, --help`: show command help.

## Examples

```bash
isamtools update sample.hg38.isam --name hg38
isamtools update sample.hg38.isam --store-path /new/sample.chm13.bam
```
