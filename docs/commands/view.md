# isamtools view

Reconstruct SAM/BAM from an alignment-view `.isam`.

```bash
isamtools view [options] <view.isam> [region]
```

## Options

- `-b`: output BAM.
- `-O, --output-fmt sam|bam|cram`: choose output format.
- `-h, --header`: include the header.
- `-H, --header-only`: print only the header.
- `--no-header`: suppress the header.
- `--store PATH`: override canonical store path.
- `--store-index PATH`: override store index path.
- `-T, --reference PATH`: reference FASTA for CRAM payload decoding.
- `-@, --threads N`: htslib I/O threads.
- `-L, --targets-file PATH`: BED file of target intervals.
- `--order coordinate|name`: full-scan order.
- `--subset NAME`: emit only reads in a view subset.
- `-f INT`: require all FLAG bits.
- `-F INT`: exclude records with any FLAG bits.
- `--rf INT`: require at least one FLAG bit.
- `-G INT`: exclude records with all FLAG bits.
- `-q INT`: minimum MAPQ.
- `-c, --count`: count records instead of emitting them.
- `-o, --output PATH`: output path. Defaults to stdout.

LMDB alignment views support full scans, name-ordered scans, filters, counts,
BED filtering, and region queries after `isamtools sort --interval` or
`--all`. Named subsets are implemented for LMDB views.

## Hard Clipping

The view CIGAR is authoritative. If hard clips make the view query length
shorter than the canonical payload, `view` trims canonical SEQ/QUAL only when
the CIGAR hard clips exactly account for the difference. Other length
mismatches are errors.

## Examples

```bash
isamtools view sample.hg38.isam chr1:1000000-2000000
isamtools view -O bam -o out.bam sample.hg38.isam
isamtools view -c -L targets.bed sample.hg38.isam
isamtools view -c -F 2304 -q 20 sample.hg38.isam
isamtools view --subset hifi sample.hg38.isam chr1:1-100000
```
