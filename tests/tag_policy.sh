#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-tag-policy.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF	RX:Z:canonA	MN:i:8
qB	0	chrS	20	60	8M	*	0	0	TTTTCCCC	FFFFFFFF	RX:Z:canonB
qC	0	chrS	40	60	8M	*	0	0	AAAACCCC	FFFFFFFF	RX:Z:canonC
SAM

cat > "$WORK/aln.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qA	0	chrV	101	30	8M	*	0	0	NNNNNNNN	!!!!!!!!	RX:Z:sourceA	HP:i:1	xx:Z:inline
qA	2048	chrV	150	10	4M4S	*	0	0	NNNNNNNN	!!!!!!!!	HP:i:2
qB	0	chrV	200	25	8M	*	0	0	NNNNNNNN	!!!!!!!!
qC	0	chrV	250	25	8M	*	0	0	NNNNNNNN	!!!!!!!!	RX:Z:canonC
SAM

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

"$TOOL" index -o "$WORK/store.isam" "$WORK/store.bam" >/dev/null 2>"$WORK/index.err"
grep -q 'read-level tags observed and reconstructible from canonical BAM/CRAM: .*RX' "$WORK/index.err"

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.isam" \
    --aln "$WORK/view.isam" 2>"$WORK/import.err"
grep -q 'stored-once tag conflicts found' "$WORK/import.err"
grep -q 'HP: 1 reads with conflicts' "$WORK/import.err"

qa=$("$TOOL" view --no-header "$WORK/view.isam" | awk '$1=="qA" && $2==0 { print; exit }')
qb=$("$TOOL" view --no-header "$WORK/view.isam" | awk '$1=="qB" { print; exit }')
qc=$("$TOOL" view --no-header "$WORK/view.isam" | awk '$1=="qC" { print; exit }')
grep -q $'RX:Z:sourceA' <<<"$qa"
grep -q $'HP:i:1' <<<"$qa"
grep -q $'xx:Z:inline' <<<"$qa"
grep -q $'RX:Z:canonB' <<<"$qb"
grep -q $'RX:Z:canonC' <<<"$qc"

stats=$("$TOOL" store_stats "$WORK/view.isam")
grep -q $'^RX\t1\t' <<<"$stats"
grep -q $'^HP\t1\t' <<<"$stats"
grep -q $'^xx\t1\t' <<<"$stats"

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.isam" \
    --drop-source-tags RX \
    --aln "$WORK/source_drop.isam" >/dev/null 2>"$WORK/source_drop.err"
qa_source_drop=$("$TOOL" view --no-header "$WORK/source_drop.isam" | awk '$1=="qA" && $2==0 { print; exit }')
grep -q $'RX:Z:canonA' <<<"$qa_source_drop"
! grep -q $'RX:Z:sourceA' <<<"$qa_source_drop"

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.isam" \
    --drop-tags RX \
    --aln "$WORK/drop.isam" >/dev/null 2>"$WORK/drop.err"
qa_drop=$("$TOOL" view --no-header "$WORK/drop.isam" | awk '$1=="qA" && $2==0 { print; exit }')
qb_drop=$("$TOOL" view --no-header "$WORK/drop.isam" | awk '$1=="qB" { print; exit }')
! grep -q $'RX:Z:' <<<"$qa_drop"
! grep -q $'RX:Z:' <<<"$qb_drop"

echo "tag policy OK"
