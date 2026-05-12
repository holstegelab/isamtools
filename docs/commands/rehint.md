# isamtools rehint

Rewrite an alignment view after its canonical store index has been rebuilt or
moved.

```bash
isamtools rehint [options] <view.isam>
```

## Options

- `--store PATH`: update/use canonical store path metadata.
- `--store-index PATH`: use this store index instead of the stored path.
- `-h, --help`: show command help.

## Behavior

`rehint` resolves every `read_dict.qname` against the target store index and
rewrites all view-side `store_read_id` columns, including alignment records,
view read tags, and subsets.

Every QNAME must resolve. Missing QNAMEs are errors.

## Examples

```bash
isamtools rehint sample.hg38.isam

isamtools rehint sample.hg38.isam \
  --store /new/sample.chm13.bam \
  --store-index /new/sample.chm13.isam
```
