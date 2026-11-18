# isamtools ref-index

Prepare the canonical genome reference as verified raw M5 files for fast
CRAM payload decoding. This is separate from `isamtools index`, which builds
the canonical read store's `.ssam` index.

```bash
isamtools ref-index --cram canon.cram -T reference.fa -o reference.m5
isamtools view --reference-path reference.m5 \
  --write-index -1 -o phased.bam phased_hifi_only.isam
```

`--cram` supplies canonical contig names, lengths and M5 checksums. `-T`
supplies the corresponding FASTA. A missing FASTA `.fai` is built automatically.
The output directory must be new; existing directories and symlinks are
refused. Its parent must exist, and its name must contain neither `:` nor `%`.

Preparation reads the FASTA in bounded chunks, uppercases the bases and
checks every contig's length and MD5 against the canonical header before
publishing that contig's file. Files contain sequence only, with no headers or
line breaks, and preserve literal M5 filename case. Repeated literal checksums
share one file. On failure, only this preparation's files and directory are
removed. Input CRAM and FASTA contents are unchanged.

To prepare and convert in one command:

```bash
isamtools view -T reference.fa --prepare-reference reference.m5 \
  --write-index -1 -o phased.bam phased_hifi_only.isam
```

`--prepare-reference DIR` requires a canonical CRAM, `-T`, and SAM/BAM record
output. It cannot be combined with `--reference-path`, count/header-only
queries or CRAM output. It leaves the prepared directory available for reuse;
subsequent conversions should use `--reference-path DIR`.

Ordinary FASTA coordinate decoding prints a brief performance hint. It keeps
using the requested FASTA and does not silently create a multi-gigabyte cache.
Count and header-only requests do not print this hint or prepare references.

Preparation time and working space depend on the genome size and storage.
Reuse the immutable directory across samples using the same canonical decoding
reference, including when exporting views aligned to a different genome. The
downstream caller still needs its selected alignment reference.
See [performance tuning](../performance.md) for input placement, cache/worker
sizing and profiling. `--write-index` creates the BAI during writing; it does
not require a separate `samtools index` command.
