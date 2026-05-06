# isamtools index

Build or rebuild a canonical BAM/CRAM store index.

```bash
isamtools index [options] <store.bam|store.cram>
```

## Options

- `--name NAME`: set a logical alias for the store.
- `--backend sqlite|lmdb`: choose the store-index backend. Default is `sqlite`.
- `-o, --output PATH`: write the index to `PATH` instead of `<store>.isam`.
- `-h, --help`: show command help.

## Behavior

`index` scans primary canonical records, assigns `canonical_read.read_id`, stores
QNAME and coordinate metadata, and builds QNAME/coordinate indexes.

When rebuilding an existing store index, store-level read tags and subsets are
preserved by QNAME and rewritten onto the new `read_id` values. If an annotated
QNAME is absent from the rebuilt store, the command fails.

LMDB store indexes are file-backed `.isam` databases and require an
LMDB-enabled build. LMDB creates a sibling lock file named `<path>-lock` while
the environment is in use.

LMDB output paths must not already exist.

## Examples

```bash
isamtools index sample.chm13.bam
isamtools index --name chm13 sample.chm13.cram
isamtools index -o sample.store.isam sample.chm13.bam
isamtools index --backend lmdb -o sample.chm13.bam.lmdb.isam sample.chm13.bam
```
