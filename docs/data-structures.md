# Data Structures And Command Paths

This document describes the LMDB-backed data structures used by `isamtools`
and the execution path taken by the main commands:

- `index`
- `import`
- `view` of a region
- `annotate`
- `rehint`

For the proposed tag ownership, fallback, and comparison model, see
`docs/tag-policy-design.md`.

The project is now LMDB-only. A `.isam` path is a file-backed LMDB environment
opened with `MDB_NOSUBDIR`, and LMDB creates the sibling lock file
`<path>-lock` while the environment exists.

## Design Summary

`isamtools` splits alignment storage into two layers:

- A canonical store index for read identity and canonical payload location.
- An alignment view for alignment-specific fields and secondary indexes.

The canonical store is keyed by one integer read identity, `read_id`.
Alignment views reuse that identity as `store_read_id`, so reconstructed SAM/BAM
records can be rebuilt by joining view state with canonical payload fetched from
the original BAM/CRAM.

The broad rule is:

- Store `.isam`: owns canonical read identity.
- View `.isam`: owns per-alignment state and optional view-local read annotations/subsets.
- `view`: reconstructs by reading the view, then fetching canonical payload from
  the store by `store_read_id`.

## Core In-Memory Structures

The public LMDB backend state lives in `lmdb_backend_t`.

Important fields:

- `env`: LMDB environment handle.
- `write_txn`: active write transaction, when a command batches writes.
- `meta`: metadata key/value DBI.
- Store DBIs: `canonical_by_id`, `canonical_by_qname`, plus optional
  store-side annotation tables such as `read_tag`, `subset_by_name`,
  `subset_read`, and `annotation_meta`.
- View DBIs: `read_by_id`, `read_by_qname`, `aln_by_id`,
  `aln_import_by_order`, `aln_by_coord`, `aln_unmapped`, `aln_by_read`,
  `aln_checkpoint`, `subset_by_name`, `subset_read`, `read_tag`,
  `annotation_meta`.

Two encoded record types are central:

- `lmdb_store_record_t`
  - `read_id`, `qname`, `qname_hash`
  - canonical placement hints: `tid`, `pos0`, `end0`
  - flags: `flag`, `is_unmapped`
  - BAM location hint: `bam_voffset`, `has_voffset`
- `lmdb_view_record_t`
  - `aln_id`, `store_read_id`
  - alignment state: `flag`, `tid`, `pos0`, `end0`, `mapq`
  - alignment payload fragments: `n_cigar`, `cigar`, `mate_tid`,
    `mate_pos0`, `tlen`, `aux`, `aux_len`
  - stable import order: `record_order`
  - unmapped flag: `is_unmapped`

Read paths use `lmdb_view_iter_t`, which is a transaction-scoped iterator over
one logical access pattern. Iterator kinds are:

- `LMDB_VIEW_ITER_IMPORT`: import order scan.
- `LMDB_VIEW_ITER_COORD`: coordinate scan, then unmapped tail.
- `LMDB_VIEW_ITER_NAME`: qname-ordered full scan.
- `LMDB_VIEW_ITER_REGION`: region scan using checkpoint-derived safe start.
- `LMDB_VIEW_ITER_READ_ID`: one read and all of its alignments.
- `LMDB_VIEW_ITER_SUBSET`: all alignments whose read is in a named subset.

Two higher-level wrappers sit on top of the backend:

- `store_t`: opens the canonical BAM/CRAM plus its LMDB store index.
- `view_t`: opens the alignment-view LMDB plus header and store linkage.

`view` and the C API also use `payload_cache_t`, an in-memory cache of rebuilt
canonical BAM records keyed by `read_id`, to avoid repeated BAM/CRAM fetches for
multiple alignments from the same read.

## On-Disk LMDB Layout

### Metadata DBI

Every `.isam` has a `meta` DBI containing string key/value metadata such as:

- `isam_type`: `store` or `alignment`
- `isam_version`
- `backend`
- `store_path`, `store_index_path`
- `indexed_name`, `indexed_coord`, `indexed_interval`
- `store_header_hash`, `store_sq_hash`
- `view_max_ref_span`, `view_max_query_len`
- `interval_strategy`, `checkpoint_step`
- `subset_next_id`, `annotation_next_id`
- operational stamps such as `created_at` and `rehinted_at`

### Store DBIs

Canonical store indexes use:

- `canonical_by_id`
  - key: big-endian `read_id`
  - value: encoded `lmdb_store_record_t`
- `canonical_by_qname`
  - key: `(qname_hash, qname)`
  - value: 8-byte `read_id`

Store `.isam` files may also carry post-index read-level annotations:

- `read_tag`
  - key: `(read_id, tag)`
  - value: encoded tag type/value/source
- `subset_by_name`, `subset_read`, `annotation_meta`
  - used when `annotate` targets the store directly

This gives two access paths:

- qname lookup during `import`, `annotate`, and `rehint`
- direct read-id lookup during reconstruction

### View DBIs

Alignment views use:

- `read_by_id`
  - key: `store_read_id`
  - value: encoded `(store_read_id, qname, qname_hash)`
- `read_by_qname`
  - key: `(qname_hash, qname)`
  - value: 8-byte `store_read_id`
- `aln_by_id`
  - key: `aln_id`
  - value: encoded `lmdb_view_record_t`
- `aln_import_by_order`
  - key: `(record_order, aln_id)`
  - value: 8-byte `aln_id`
- `aln_by_coord`
  - key: `(tid, pos0, record_order, aln_id)`
  - value: 8-byte `aln_id`
- `aln_unmapped`
  - key: `(record_order, aln_id)`
  - value: 8-byte `aln_id`
- `aln_by_read`
  - key: `(store_read_id, record_order, aln_id)`
  - value: 8-byte `aln_id`
- `aln_checkpoint`
  - key: `(tid, checkpoint_pos)`
  - value: earliest safe coordinate scan start for that checkpoint
- `subset_by_name`
  - key: subset name
  - value: encoded subset metadata including `subset_id`
- `subset_read`
  - key: `(subset_id, read_id)`
  - value: empty
- `read_tag`
  - key: `(read_id, tag)`
  - value: encoded tag type/value/source
- `annotation_meta`
  - key: `annotation_id`
  - value: encoded annotation provenance

### Ordering And Scan Tables

Several fields and DBIs exist only to make scan order deterministic and cheap.

- `record_order`
  - assigned once during `import` as the zero-based kept-record ordinal
  - preserves a stable full-scan order identical to input encounter order
  - acts as a deterministic tie-breaker when multiple alignments share the same
    coordinate or the same `store_read_id`
- `aln_import_by_order`
  - is the append-oriented full-scan table used before coordinate indexes exist
  - gives a stable fallback for sequential reads and for rebuilding secondary
    indexes later in `sort`
- `aln_unmapped`
  - holds unmapped alignments separately from `aln_by_coord`
  - is used by coordinate full scans so output order is mapped records in
    coordinate order followed by unmapped records in import order

Why `record_order` is present in other keys:

- in `aln_by_coord`, it breaks ties for alignments with the same `(tid, pos0)`
  without having to inspect payloads
- in `aln_by_read`, it preserves the original encounter order of all
  alignments for one read, including primary, secondary, and supplementary
  records
- in both DBIs, it helps keep output deterministic across rebuilds of the
  secondary indexes

Two important design choices shape the rest of the system:

- `aln_by_id` owns the full alignment payload.
- Secondary DBIs mostly store compact 8-byte references to `aln_id`.

That keeps the ordered indexes small and reduces rewrite cost.

## Record Encoding Rules

The LMDB keys and values are manually encoded rather than storing raw C
structs.

Key properties:

- Integers are written big-endian, so LMDB lexical ordering matches numeric
  ordering.
- Strings and blobs are length-prefixed.
- Values carry a record version byte.
- Unknown or malformed encodings are fatal errors.

This matters for ordered scans. For example:

- `canonical_by_qname` sorts by `qname_hash` then `qname`.
- `aln_by_coord` sorts first by reference, then coordinate, then stable import
  order.
- `aln_by_read` groups every alignment belonging to the same read together.

## Command Paths

### `index`

Purpose: build the canonical store index for a BAM or CRAM.

Primary input:

- canonical BAM/CRAM file

Primary output:

- store `.isam` LMDB

Execution path:

1. Open the canonical file with htslib.
2. Read the SAM header.
3. Create the LMDB environment and begin a write transaction.
4. Write store-level metadata into `meta`.
5. Derive and store two header hashes:
   - `store_header_hash` from raw header text
   - `store_sq_hash` from the `@SQ` projection
6. Stream the canonical file record by record.
7. Skip secondary and supplementary records so one canonical identity is kept
   per read name.
8. Assign `read_id = n_records + 1` in encounter order.
9. Determine the effective read-level tag list from built-in tag lists, header
  hints, and any `--read-level-tags` overrides.
10. Encode and insert the record into:
   - `canonical_by_id`
   - `canonical_by_qname`
11. Store the effective read-level tag list in `meta`.
12. Keep canonical aux tags in the original BAM/CRAM payload rather than
  duplicating them into store LMDB.
13. Track `store_max_read_len` and, for BAM, record `bam_voffset` as a payload
    hint.
14. Periodically commit and reopen the write transaction for large inputs.
15. Write final index metadata such as `indexed_name` and `indexed_coord`.

Why the two store indexes exist:

- `canonical_by_qname` is the identity resolver used by `import`, `annotate`,
  and `rehint`.
- `canonical_by_id` is the direct lookup used later when reconstructing reads.

Canonical-store consequence:

- `store_read_id` always points at the single canonical store record for the
  read, not at a supplementary alignment
- this is intentional because supplementary records may be hard-clipped or
  otherwise unsuitable as the canonical payload source
- the current implementation prefers a primary alignment and otherwise keeps a
  non-hard-clipped supplementary alignment as the canonical payload source
- canonical BAM/CRAM read-level tags stay in the canonical payload; store
  `read_tag` is used only for later store-side annotations

### `import`

Purpose: build an alignment view that points back to a canonical store.

Primary inputs:

- alignment SAM/BAM/CRAM
- canonical store BAM/CRAM
- canonical store `.isam`

Primary output:

- alignment-view `.isam` LMDB

Execution path:

1. Open the canonical store via `store_open`.
2. Verify that the store has a qname index.
3. Open the input alignment stream and read its header.
4. Create the output view LMDB.
5. Write view metadata:
   - type/version/backend
   - store path and store index path
   - copied store header hashes and max read length
   - index availability flags, initially disabled
6. Store the alignment-view header text in `meta`.
7. Begin a batched LMDB write transaction.
8. For each input alignment:
   - extract `qname`
   - look up the canonical read in `canonical_by_qname`
   - fail or skip if absent, depending on `--allow-missing`
   - classify aux tags using the stored read-level tag list, the input header,
     and any import-time overrides
   - promote read-level tags into `read_tag`
   - keep alignment-local and unknown tags in the inline aux blob
   - write a read identity entry into `read_by_id` and `read_by_qname`
   - encode the alignment into `aln_by_id`
   - add import-order and read fanout entries to:
     - `aln_import_by_order`
     - `aln_by_read`
9. Track `view_max_ref_span` and `view_max_query_len` for later region logic.
10. Commit in batches for large imports.
11. Optionally invoke `sort --all`.

What `import` intentionally does not build:

- `aln_by_coord`
- `aln_unmapped`
- `aln_checkpoint`

Those are built later by `sort` so import stays append-oriented.

Tag consequence during import:

- read-level SAM tags live in `read_tag`
- alignment-local SAM tags stay inline in each `lmdb_view_record_t` as the raw
  aux blob `aux` plus `aux_len`
- unknown tags default to alignment-level and are reported as warnings
- this keeps the immutable alignment-view model while still letting
  reconstruction merge read-level tags back into output

### `view` of a Region

Purpose: reconstruct records overlapping one genomic interval.

Primary inputs:

- alignment-view `.isam`
- linked canonical store BAM/CRAM plus store `.isam`
- optional BED file or subset filter layered on top

Primary output:

- reconstructed SAM/BAM/CRAM stream or a count

Execution path for a plain region query:

1. Open the view with `view_open_read`.
2. Resolve and open the linked canonical store with `view_link_store`.
3. Parse the region against the view header.
4. Require `indexed_interval=1`. Region queries assume `sort --interval` or
   `sort --all` has already run.
5. Open `LMDB_VIEW_ITER_REGION`.
6. Inside iterator setup:
   - open a read transaction
   - open a cursor on `aln_by_coord`
   - compute `scan_start`
7. `scan_start` comes from one of two strategies:
   - `max_span`: scan backward by `view_max_ref_span`
   - `checkpoint`: read `aln_checkpoint` to find the earliest still-active
     alignment start near the query boundary
8. Seek `aln_by_coord` to `(tid, scan_start, 0, 0)`.
9. Walk forward in coordinate order until `tid` changes or `pos0 >= region_end`.
10. For each candidate alignment:
    - dereference `aln_id` into `aln_by_id`
    - exact-filter on overlap: `pos0 < end && end0 > beg`
11. For each passing record, reconstruct the full BAM record:
    - resolve `store_read_id` in `canonical_by_id`
    - fetch canonical payload from BAM/CRAM, preferably via cached voffset
    - rebuild a fresh `bam1_t` with view-side alignment state
    - append the stored view-side aux blob, which includes the
      alignment-specific SAM tags imported from the source alignment file
    - merge view `read_tag` entries for the same `read_id`
    - merge store `read_tag` entries for the same `read_id` when not already
      present in the record; these are store-side annotations rather than a
      duplicate copy of canonical BAM/CRAM tags
12. Apply optional flag/MAPQ/BED filters and emit or count.

Why this is structured as a two-stage lookup:

- `aln_by_coord` is optimized for ordered search.
- `aln_by_id` owns the full record payload.
- Exact interval filtering happens after dereference so long-read overlap is
  correct even when scanning from an earlier safe start.

Important current behavior:

- `view` rebuilds alignment-local tags from the inline `aln_by_id` aux blob
- `read_tag` entries from the view are merged back into the emitted SAM/BAM
  record during reconstruction
- store `read_tag` entries are also merged back, but those are only later
  store-side annotations; canonical BAM/CRAM tags already came from the fetched
  payload
- view-level read tags take precedence over store-level read tags when both are
  present for the same `(read_id, tag)` pair

Special cases layered onto the same base path:

- `--subset` plus region: iterate the subset and then filter overlaps in memory.
- `-L/--targets-file`: run multiple region iterators and deduplicate by `aln_id`.

### `annotate`

Purpose: attach read-level tags or add reads to named subsets.

Primary inputs:

- target store `.isam` or view `.isam`
- TSV tag file or one-qname-per-line subset file

Primary output:

- updated target `.isam`

Two modes exist.

#### Tag mode

1. Open the target LMDB.
2. Inspect `isam_type` to decide whether qname resolution goes through:
   - `canonical_by_qname` for a store
   - `read_by_qname` for a view
3. Read the input TSV line by line.
4. For each qname, resolve the integer read identity.
5. For each `TAG:TYPE:VALUE` field on the line, upsert one `read_tag` entry
   keyed by `(read_id, tag)`.
6. Append an `annotation_meta` record describing provenance.

#### Subset mode

1. Open the target LMDB.
2. Resolve or create the subset using `subset_by_name`.
3. If creating a new subset, allocate `subset_id` from `meta.subset_next_id`.
4. Read qnames from the subset file.
5. Resolve each qname to a read identity.
6. Insert `(subset_id, read_id)` into `subset_read` with `MDB_NOOVERWRITE`.
7. Append an `annotation_meta` record.

Why annotations are read-level, not alignment-level:

- The identity model is anchored on canonical reads.
- Subsets and tags follow the read across supplementary/secondary alignments.
- `rehint` can then remap annotation keys alongside `store_read_id` changes.

Store-side annotation status:

- annotating a store `.isam` is supported
- `cmd_annotate` resolves qnames through `canonical_by_qname` when
  `isam_type=store`
- this makes store-side subset and read-tag writes cheap
- reconstruction applies read tags from the alignment view first and then from
  the store
- `view --subset` searches the alignment view first and then falls back to the
  store subset table

### `rehint`

Purpose: rewrite stale view-side `store_read_id` values after the canonical
store index has been rebuilt and read IDs have changed.

Primary inputs:

- alignment-view `.isam`
- target store `.isam`
- optionally updated `--store` and `--store-index` metadata

Primary output:

- same view `.isam`, with rewritten read identities and dependent indexes

Execution path:

1. Open the view LMDB read-write.
2. Verify `isam_type=alignment`.
3. Resolve the target store path and store index path from either CLI overrides
   or existing view metadata.
4. Open the target store LMDB read-only and verify `isam_type=store`.
5. Start a read transaction on the store and a write transaction on the view.
6. Scan `read_by_id` in the view.
7. For each entry:
   - decode `(old_id, qname, qname_hash)`
   - resolve the qname in the target store’s `canonical_by_qname`
   - build a sorted mapping `old_id -> new_id`
8. Scan and stage dependent read-keyed structures before destructive changes:
   - `subset_read`
   - `read_tag`
9. Drop and rebuild the read-keyed secondary DBIs:
   - `read_by_id`
   - `read_by_qname`
   - `aln_by_read`
   - `subset_read`
   - `read_tag`
10. Rewrite every `aln_by_id` record in place so `store_read_id` points at the
    new store identity.
11. Recreate `aln_by_read` from the rewritten records.
12. Reinsert remapped subset membership and tags using the new read IDs.
13. Update view metadata:
    - `store_path`
    - `store_index_path`
    - `store_header_hash`
    - `store_sq_hash`
    - `rehinted_at`
    - `rehinted_read_count`
14. Commit atomically.

What `rehint` does not need to rebuild:

- `aln_by_coord`
- `aln_unmapped`
- `aln_import_by_order`
- `aln_checkpoint`

Those DBIs are keyed by alignment position or stable import/alignment identity,
not by `store_read_id`, so they remain valid.

## Why Region Queries Need `sort --interval`

`import` only writes:

- the alignment payload table (`aln_by_id`)
- import-order scan (`aln_import_by_order`)
- read fanout (`aln_by_read`)

`sort --coord` adds:

- `aln_by_coord`
- `aln_unmapped`

`sort --interval` adds:

- `aln_checkpoint`
- region metadata such as `interval_strategy` and `checkpoint_step`

So the minimum path for an efficient region view is:

1. `index`
2. `import`
3. `sort --interval` or `sort --all`
4. `view <region>`

## Practical Read/Write Ownership

If you need to know where a command mutates data, this is the short version:

- `index`
  - writes: store `meta`, `canonical_by_id`, `canonical_by_qname`
- `import`
  - writes: view `meta`, `read_by_id`, `read_by_qname`, `aln_by_id`,
    `aln_import_by_order`, `aln_by_read`
- `sort --coord`
  - writes: `aln_by_coord`, `aln_unmapped`
- `sort --interval`
  - writes: `aln_checkpoint`
- `annotate --from`
  - writes: `read_tag`, `annotation_meta`
- `annotate --subset`
  - writes: `subset_by_name`, `subset_read`, `annotation_meta`
- `rehint`
  - rewrites: `read_by_id`, `read_by_qname`, `aln_by_id`, `aln_by_read`,
    `subset_read`, `read_tag`, and related metadata

## Typical End-To-End Flow

For the common case of importing one alternate alignment and querying regions:

1. `index store.bam`
2. `import aln.bam --store store.bam --aln aln.isam`
3. `sort --all aln.isam`
4. `annotate` if needed for tags or named subsets
5. `view aln.isam chr:start-end`
6. If the store index is rebuilt later, run `rehint aln.isam`

That sequence matches the current storage model: append first, build ordered
indexes second, annotate at the read level, and repair identity links only when
the canonical store’s integer IDs change.