# isamtools query

Build a QNAME selection and emit TSV fields from one or more named sources.

```bash
isamtools query \
  --source hg38:reads.hg38.ssam \
  --source chm13:reads.chm13.isam \
  --select hg38:chrX:1000000-1100000 \
  --columns qname,hg38.region,chm13.region,chm13.in_subset:hifi
```

## Sources

- `--source ALIAS:PATH`: add a named source.
- `PATH` may be a store `.ssam`, an alignment-view `.isam`, or a store
  BAM/CRAM/SAM with the matching `.ssam`.
- The implicit alias `ssam` resolves to the only store source, and `isam`
  resolves to the only alignment-view source. If either is ambiguous, use an
  explicit alias.
- Store sources report the canonical store alignment only. Alignment views may
  report multiple alignments per read; list-valued fields are joined with `;`.

## Read Selection

- `--select ALIAS:REGION`, `--region ALIAS:REGION`: select reads overlapping a
  source-specific region.
- `--subset ALIAS:NAME`: select reads from a named subset in that source. For
  alignment views, store-level subsets are also accepted.
- `--qnames FILE`: select one QNAME per line.
- `--join union|inner`: combine selectors by union or intersection. The default
  is `union`.

Region selection on alignment views requires `isamtools sort --interval` or
`--all`. Region selection on store sources uses the canonical BAM/CRAM index.

## Columns

Columns are comma-separated. `qname` is global. Source columns use
`ALIAS.FIELD`.

Supported source fields:

- `read_id`
- `n` or `n_aln`
- `rname`
- `pos`
- `start0`
- `end0`
- `region`
- `flag`
- `mapq`
- `strand`
- `cigar`
- `is_mapped`
- `in_subset:NAME`

Coordinates in `pos` are 1-based SAM-style starts. `start0` and `end0` are
0-based half-open coordinates. `region` is printed as `rname:pos-end`, with
1-based inclusive-looking endpoints matching the aligned span.

## Examples

```bash
isamtools query \
  --source hg38:reads.hg38.ssam \
  --source chm13:reads.chm13.isam \
  --select hg38:chrX:1000000-1100000 \
  --columns qname,hg38.region,hg38.mapq,chm13.region,chm13.mapq

isamtools query \
  --source hg38:reads.hg38.ssam \
  --source chm13:reads.chm13.isam \
  --select hg38:chrX:1000000-1100000 \
  --select chm13:chrX:2000000-2100000 \
  --join inner \
  --columns qname,hg38.region,chm13.region
```
