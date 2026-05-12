#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-lmdb-aln.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
qB	0	chrS	20	50	8M	*	0	0	TTTTCCCC	FFFFFFFF
qC	4	*	0	0	*	*	0	0	GGGGTTTT	FFFFFFFF
qD	0	chrS	40	50	8M	*	0	0	CCCCAAAA	FFFFFFFF
SAM

cat > "$WORK/aln.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qA	0	chrV	101	30	8M	*	0	0	NNNNNNNN	!!!!!!!!	NM:i:0
qB	16	chrV	200	20	8M	*	0	0	NNNNNNNN	!!!!!!!!	NM:i:1
qA	2048	chrV	150	10	4M4S	*	0	0	NNNNNNNN	!!!!!!!!	SA:Z:chrV,101,+,8M,30,0;
qC	4	*	0	0	*	*	0	0	NNNNNNNN	!!!!!!!!
qD	0	chrV	300	30	8M	*	0	0	NNNNNNNN	!!!!!!!!
SAM

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

"$TOOL" index -o "$WORK/store.lmdb.isam" "$WORK/store.bam"
"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.lmdb.isam" \
    --aln "$WORK/view.lmdb.isam"
test -f "$WORK/view.lmdb.isam"
test -s "$WORK/view.lmdb.isam"
test -f "$WORK/view.lmdb.isam-lock"

[[ $("$TOOL" view -c "$WORK/view.lmdb.isam") == 5 ]]
[[ $("$TOOL" view --no-header --order name "$WORK/view.lmdb.isam" | awk 'END { print NR }') == 5 ]]

"$TOOL" sort --all "$WORK/view.lmdb.isam"

[[ $("$TOOL" view -c "$WORK/view.lmdb.isam") == 5 ]]
[[ $("$TOOL" view -c -f 2048 "$WORK/view.lmdb.isam") == 1 ]]
[[ $("$TOOL" view -c -F 2048 "$WORK/view.lmdb.isam") == 4 ]]
[[ $("$TOOL" view --no-header "$WORK/view.lmdb.isam" chrV:145-155 | awk 'END { print NR }') == 1 ]]

seq_a=$("$TOOL" view --no-header "$WORK/view.lmdb.isam" \
    | awk '$1=="qA" && $2==0 { print $10; exit }')
[[ $seq_a == ACGTACGT ]]

"$TOOL" update "$WORK/view.lmdb.isam" --name lmdb-view

echo "lmdb alignment OK"