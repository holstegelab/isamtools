# isamtools index

Build or rebuild a canonical BAM/CRAM store index.

```bash
isamtools index [options] <store.bam|store.cram>
```

## Options

- `--name NAME`: set a logical alias for the store.
- `-o, --output PATH`: write the index to `PATH` instead of `<store-root>.isam`.
- `--share-tags TAGS`: add comma-separated aux tags to canonical
  reconstruction from the BAM/CRAM payload.
- `--no-share-tags TAGS`: remove comma-separated aux tags from canonical
  reconstruction.
- `--strict-tag-policy`: fail when unknown canonical tags are unresolved.
- `--verbose-tag-policy`: print the detailed tag policy summary.
- `--yes`: accept tag policy warnings.
- `-h, --help`: show command help.

## Behavior

`index` scans primary canonical records, assigns `canonical_read.read_id`, stores
QNAME and coordinate metadata, and builds QNAME/coordinate indexes.

The default output path replaces a trailing `.bam`, `.cram`, or `.sam` with
`.isam`. For example, `sample.chm13.bam` becomes `sample.chm13.isam`.

Primary canonical records must not contain hard clipping. The canonical
BAM/CRAM is the source of full SEQ/QUAL, so hard-clipped reads belong in
alignment views, not in the store.

Store indexes are file-backed LMDB `.isam` databases. LMDB creates a sibling
lock file named `<path>-lock` while the environment is in use.

LMDB output paths must not already exist.

The index records observed canonical reconstructible tags in store metadata.
Those aux values are not copied into the `.isam`; they are read from the
canonical BAM/CRAM payload during `view` when an imported alignment file does
not provide its own value. Unknown canonical tags are reported but are not
reconstructed unless selected with `--share-tags`.

`index` prints the intended tag-policy behavior before scanning, then prints an
observed tag summary at the end. The observed summary only lists tags that were
actually seen.

## Examples

```bash
isamtools index sample.chm13.bam
isamtools index --name chm13 sample.chm13.cram
isamtools index -o sample.store.isam sample.chm13.bam
isamtools index --share-tags RX,MM,ML --no-share-tags MN sample.chm13.bam
```
