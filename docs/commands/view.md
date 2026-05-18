# isamtools view

Reconstruct SAM/BAM from an alignment-view `.isam`.

```bash
isamtools view [options] <view.isam> [region]
```

## Options

- `-b, --bam`: output BAM.
- `-C, --cram`: output CRAM; requires `-T`.
- `-1, --fast`: fast BAM compression and default output to BAM.
- `-u, --uncompressed`: uncompressed BAM output.
- `-O, --output-fmt sam|bam|cram[,level=N|fast|uncompressed]`: choose output format.
- `--output-fmt-option level=N|fast|uncompressed`: output compression option.
- `-h, --header, --with-header`: include the header in SAM output.
- `-H, --header-only`: print only the header.
- `--no-header`: suppress SAM header output. This is the default for SAM, matching `samtools view`.
- `--store PATH`: override canonical store path.
- `--store-index PATH`: override store index path.
- `-T, --reference PATH`: reference FASTA for CRAM payload decoding and CRAM output.
- `-@, --threads N`: htslib I/O threads.
- `-L, --targets-file PATH`, `--regions-file PATH`: BED file of target intervals.
- `-M, --use-index`: accepted for `samtools view` compatibility; indexed paths are used automatically.
- `-P, --fetch-pairs`: for region or BED selection, expand selected read names to all alignments for those reads.
- `--order coordinate|name`: full-scan order.
- `--subset NAME`: emit only reads in a named subset, searching the alignment view first and then the canonical store.
- `-f INT`: require all FLAG bits.
- `-F INT`: exclude records with any FLAG bits.
- `--rf INT`: require at least one FLAG bit.
- `-G INT`: exclude records with all FLAG bits.
- `-q INT`: minimum MAPQ.
- `-m INT`: minimum query bases consumed by the CIGAR.
- `-r STR`: require read group `RG:Z:STR`.
- `-R PATH`: require read group listed in `PATH`.
- `-N PATH`: require QNAME listed in `PATH`.
- `-l STR`: require library `STR`, resolved through the output header `@RG ID -> LB`.
- `-d TAG[:VALUE]`: require an aux tag, optionally with a matching value.
- `-D TAG:PATH`: require an aux tag whose value appears in `PATH`.
- `--subsample F`: keep fraction `F` of reads, from `0` to `1`.
- `--subsample-seed INT`: seed used by `--subsample`.
- `-s INT.FRAC`: samtools-compatible shorthand for `--subsample 0.FRAC --subsample-seed INT`.
- `-x TAGS`, `--remove-tag TAGS`: strip comma-separated aux tags from output.
  `-x ^TAGS` preserves only those tags, matching samtools shorthand.
- `--keep-tag TAGS`: preserve only comma-separated aux tags in output.
- `-p`, `--unmap`: emit records not selected by filters/regions with the
  `UNMAP` flag set instead of dropping them.
- `-c, --count`: count records instead of emitting them.
- `-o, --output PATH`: output path. Defaults to stdout.
- `-S`, `--no-PG`, and `-?` are accepted for `samtools view` compatibility.

LMDB alignment views support full scans, name-ordered scans, filters, counts,
BED filtering, one or more positional regions, and region queries after
`isamtools sort --interval` or `--all`. Named subsets are searched in the
alignment view first and then in the linked canonical store.

## Hard Clipping

Primary canonical reads in the store must not be hard-clipped. Hard clipping is
allowed in alignment views. The view CIGAR is authoritative; if hard clips make
the view query length shorter than the canonical payload, `view` trims
canonical SEQ/QUAL only when the hard clips exactly account for the difference.
Other length mismatches are errors.

## Examples

```bash
isamtools view sample.hg38.isam chr1:1000000-2000000
isamtools view -h sample.hg38.isam chr1:1000000-2000000
isamtools view -b -1 -o out.bam sample.hg38.isam
isamtools view -C -T ref.fa -o out.cram sample.hg38.isam
isamtools view -c -L targets.bed sample.hg38.isam
isamtools view -c -N qnames.txt -r rg1 sample.hg38.isam chr1:1-100000 chr2:1-100000
isamtools view -P sample.hg38.isam chr1:1000000-2000000
isamtools view -c -F 2304 -q 20 sample.hg38.isam
isamtools view --subsample 0.1 --subsample-seed 42 sample.hg38.isam
isamtools view -x MM,ML sample.hg38.isam
isamtools view -p sample.hg38.isam chr1:1000000-2000000
isamtools view --subset hifi sample.hg38.isam chr1:1-100000
```
