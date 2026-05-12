#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-tag-audit.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF	rt:Z:readA	zz:Z:store
qB	0	chrS	20	60	8M	*	0	0	TTTTCCCC	FFFFFFFF	rt:Z:readB
SAM

cat > "$WORK/aln.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qA	0	chrV	101	30	8M	*	0	0	NNNNNNNN	!!!!!!!!	rt:Z:readA	NM:i:0	zz:Z:aln
qA	2048	chrV	150	10	4M4S	*	0	0	NNNNNNNN	!!!!!!!!	rt:Z:readA	NM:i:2	zz:Z:split
qB	0	chrV	200	25	8M	*	0	0	NNNNNNNN	!!!!!!!!	rt:Z:readB	NM:i:1
SAM

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

"$TOOL" index --share-tags rt -o "$WORK/store.isam" "$WORK/store.bam" 2>"$WORK/index.err"
grep -q 'unknown canonical aux tags are not shared: zz' "$WORK/index.err"

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.isam" \
    --aln "$WORK/view.isam" 2>"$WORK/import.err"
grep -q 'unknown aux tags default to alignment retention: zz' "$WORK/import.err"

"$TOOL" sort --all "$WORK/view.isam"

line=$("$TOOL" view --no-header "$WORK/view.isam" | awk '$1=="qA" && $2==0 { print; exit }')
grep -q $'rt:Z:readA' <<<"$line"
grep -q $'NM:i:0' <<<"$line"
grep -q $'zz:Z:aln' <<<"$line"
[[ $(grep -o 'rt:Z:readA' <<<"$line" | wc -l) -eq 1 ]]

audit=$("$TOOL" audit-tags --share-tags rt "$WORK/aln.sam")
grep -q $'^tag\tscope\tsource\talignments$' <<<"$audit"
grep -q $'^NM	alignment	builtin-alignment	3$' <<<"$audit"
grep -q $'^rt	shared	configured-shared	3$' <<<"$audit"
grep -q $'^zz	unknown	unknown	2$' <<<"$audit"

echo "tag audit OK"
