---
name: isamtools
description: >-
  Use isamtools to deduplicate multi-alignment sequencing data: keep ONE canonical
  BAM/CRAM (the only copy of SEQ/QUAL) plus lightweight `.isam` alignment views for every
  other alignment of the same reads (another reference such as hg38 vs CHM13, a re-alignment,
  a phased/filtered variant). Use when building, importing into, querying, validating,
  archiving, moving or restoring isamtools bundles; when reconstructing BAM/SAM/CRAM from a
  view; when deciding whether original BAMs can be deleted; or when debugging tag-policy,
  CRAM-reference or missing-read errors. Triggers: "isamtools", ".isam", ".ssam",
  "canonical store", "alignment view", "dedup BAMs across references", "replace BAMs with
  isam", "reconstruct BAM from isam".
---

# isamtools

`isamtools` stores read payload (SEQ, QUAL and read-level tags such as PacBio `MM`/`ML`)
**once**, in a canonical BAM or CRAM. Every other alignment of the same reads becomes an
`.isam` view (LMDB) that keeps only alignment fields (FLAG, RNAME, POS, MAPQ, CIGAR, mate
fields, alignment-level tags). `isamtools view` rebuilds full records on the fly.

A complete **bundle** is:

```
canon.bam | canon.cram (+ .bai/.crai)   canonical payload (only copy of SEQ/QUAL)
canon.ssam                              store index: QNAME -> payload locator, tag fingerprints
<name>.isam  (one per alignment)        alignment views
```

Losing the canonical file or its `.ssam` makes **every** view unreadable. Treat the
bundle as one unit when copying, archiving or deleting.

## Build

```bash
make                       # vendors + statically links htslib, LMDB, bzip2
make test                  # functional tests (needs samtools, Python 3.9+)
```

- Build on the machine (or an ABI-compatible one) where it will run; a binary built
  against a newer glibc will not start on an older system. Use a separate build directory
  per machine: `make BUILD_DIR=build-<host>`.
- After editing any header, rebuild from clean (`make BUILD_DIR=... clean`). A stale
  object compiled against an old struct layout silently reads garbage.

## Core workflow

```bash
# 1. Canonical store. Pick the alignment holding the UNION of all reads (e.g. a combined
#    HiFi+non-HiFi BAM), primary records not hard-clipped.
samtools view -O CRAM -T ref.fa -o canon.cram canonical.bam   # optional: CRAM canonical
samtools index canon.cram
isamtools index -T ref.fa canon.cram                          # -> canon.ssam

# 2. One view per additional alignment of the same reads
isamtools import --store canon.cram --store-index canon.ssam \
    --aln other.isam --strict-tag-policy --fast-verify other.bam
isamtools sort --all other.isam                               # coordinate/name/interval indexes

# 3. Make the bundle relocatable (paths resolved relative to the .isam)
isamtools update other.isam --store-path canon.cram --store-index-path canon.ssam

# 4. Use it like a BAM
isamtools view -T ref.fa other.isam chr1:1000000-2000000
isamtools view -b -T ref.fa -o other.bam other.isam           # full BAM back
```

`import` reads from `-` (stdin), so a source BAM can be streamed from a remote host without
a local copy: `ssh host cat x.bam | isamtools import --store canon.cram ... -`.

### Reference rules (CRAM canonical)

- Every view decodes payload from the **canonical** file, so `-T` is always the
  **canonical's** reference, even for views aligned to a different genome (a CHM13 view
  over an hg38 CRAM canonical still needs the hg38 FASTA).
- The FASTA must be the exact one used to write the CRAM (`@SQ M5` checksums must match).
- For fast whole-file export, prepare verified raw references once with `isamtools ref-index`
  (see `docs/commands/ref-index.md` and `docs/performance.md`).

## Before deleting the original BAMs: validate properly

`import` stores alignment fields and checks shared tags (`--fast-verify` compares
per-tag fingerprints, `--strict-tag-policy` aborts on conflicts). **It does not compare
SEQ/QUAL**: those are looked up from the canonical by QNAME at view time. Record-count
equality is therefore necessary but not sufficient. Before removing a source BAM:

1. **Counts:** `samtools idxstats`/`view -c` of the source == `isamtools view -c` of the
   view, including unmapped (`-f 4`). `import` must report `0 missing`.
2. **Exact content:** `isamtools check --view other.isam other.bam` compares every record
   (SEQ, QUAL, all tags). For the canonical, `isamtools check --store canon.cram
   --store-index canon.ssam canonical.bam` compares **primary** records only (secondary and
   supplementary are reported as skipped), so also compare the canonical CRAM itself with
   the source BAM, e.g. by streaming both through `samtools view` into a normaliser.
3. **Independent spot check** (catches bugs shared by import and check): extract the same
   regions plus a slice of unmapped reads (`'*'`) from the source with `samtools view`
   and from the bundle with `isamtools view`, normalise, and compare every field with a
   separate script. Include supplementary/hard-clipped and alt-contig regions. Add a
   negative control (mutate one base, confirm the comparison fails).
4. **Headers:** `samtools view -H` of source vs `isamtools view -H` of the view.
5. **Transport:** checksum every bundle file at the destination (sha256/Adler-32) before
   deleting anything at the source.

Delete only the exact source files a bundle was built from. Other files in the same
directory (unaligned/per-run BAMs, assembly-contig BAMs, unphased intermediates) are **not**
covered by the bundle.

### Expected, value-preserving differences

These are not losses; normalise them in comparisons:

| Difference | Cause |
|---|---|
| CIGAR `=`/`X` come back as `M` (same lengths) | htslib CRAM decoding of a CRAM canonical |
| unmapped records: MAPQ `255` comes back as `0` | htslib CRAM decoder sentinel |
| integer tag subtype changes (`C` vs `i`), same value | BAM integer storage width |
| extra `MD`/`NM` on CRAM-decoded canonical records | CRAM decoder regenerates them |
| extra `@SQ M5:/UR:` header fields | CRAM writer |
| extra policy-derived tags in reconstructed views | shared-tag fallback |

`isamtools check` already accepts exactly these; any other difference fails.

### Interpreting import warnings

- `stored-once tag conflicts found; first observed values were retained`: without
  `--strict-tag-policy` this can mean a real per-read value difference was collapsed. Use
  `--strict-tag-policy` for archival imports, or prove losslessness with `check`.
- `QNAME '...' not found in canonical store`: the canonical does not contain that read
  (wrong canonical, or only hard-clipped supplementary alignments). Do not paper over it
  with `--allow-missing` for an archive.
- `index: N reads had only hard-clipped supplementary alignments and were not indexed`:
  those reads cannot back any view; choose a canonical with full primaries.
- Systematically differing CCS scalar tags (`ec`, `rq`, `sn`, `np`, `zm`) between two
  processing runs of the same reads are stored per read automatically; reconstruction
  stays exact.

## Archiving, moving, restoring

- Run `isamtools update ... --store-path canon.cram --store-index-path canon.ssam` before
  archiving, so a restored bundle works from any directory. Views with absolute store paths
  still work after restore if you pass `--store canon.cram --store-index canon.ssam`
  (relative to the caller's working directory).
- Archive every file plus a checksum manifest and a completion marker written **last**, so a
  partially uploaded bundle is recognisable.
- `.isam`/`.ssam` are LMDB environments. Do not open the same file from several nodes
  concurrently on a shared/parallel filesystem (reader locks are not cross-node safe); copy
  the small `.isam`/`.ssam` to node-local storage per reader. The large CRAM can stay shared.
- `*-lock` files beside `.isam`/`.ssam` are LMDB lock files; they need not be archived.

## Tools that only read BAM (pysam, unmodified programs)

- Recompile htslib-based programs with the shim (`docs/htslib-compat-migration.md`), or
- materialise a temporary indexed BAM for the run: `tools/run_with_bam.py` builds and
  checks the BAM on local scratch, runs the consumer, and deletes the temporary file
  (`docs/temporary-bam.md`, tuning in `docs/performance.md`).

## Command map

| Command | Use |
|---|---|
| `index` | build `.ssam` over the canonical BAM/CRAM (`-T` for CRAM, `--share-tags`) |
| `import` | add an alignment as `.isam` (`--strict-tag-policy`, `--fast-verify`, `--check`) |
| `sort --all` | build view indexes (needed for region queries) |
| `view` | reconstruct SAM/BAM/CRAM; samtools-like filters (`-f/-F/-q/-N/-r/-L`, `'*'` = unmapped) |
| `check` | exact record comparison of a BAM/CRAM/SAM against a view or store |
| `update` | rewrite store paths/aliases (portable bundles) |
| `annotate` | add per-read tags or named QNAME subsets (`view --subset`) |
| `audit-tags` | preview how aux tags will be classified |
| `query` | multi-source per-read TSV fields |
| `store_stats` | storage breakdown of a `.ssam` |
| `rehint` | repair store references after re-indexing a store |
| `ref-index` | prepare checksum-verified raw references for fast CRAM decoding |

Per-command details: `docs/commands/<command>.md`. Layout: `docs/data-structures.md`.
