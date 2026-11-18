# isamtools update

Mutate metadata in a store `.ssam` or alignment-view `.isam`.

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

For a movable archive bundle with the canonical files beside the view:

```bash
isamtools update bundle/sample.hg38.isam \
  --store-path canon.cram --store-index-path canon.ssam
```

`update` stores these metadata values verbatim. Relative store paths are
resolved by current readers from the directory containing the `.isam` (or
`.ssam` for its own `store_path`), not from the caller's working directory.
The bundle therefore remains readable after moving it, without rewriting the
bulk data. Explicit `view --store` / `--store-index` overrides still use the
caller's working directory. Older readers need overrides or a working directory
inside the bundle.
