#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-lmdb-store.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
r1	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
r2	0	chrS	20	50	8M	*	0	0	TTTTCCCC	FFFFFFFF
r3	4	*	0	0	*	*	0	0	GGGGTTTT	FFFFFFFF
SAM

cat > "$WORK/view.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
r1	0	chrV	101	30	8M	*	0	0	NNNNNNNN	!!!!!!!!
r2	16	chrV	200	20	8M	*	0	0	NNNNNNNN	!!!!!!!!
r3	4	*	0	0	*	*	0	0	NNNNNNNN	!!!!!!!!
SAM

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

"$TOOL" index -o "$WORK/store.lmdb.ssam" "$WORK/store.bam"
test -f "$WORK/store.lmdb.ssam"
test -s "$WORK/store.lmdb.ssam"
test -f "$WORK/store.lmdb.ssam-lock"

"$TOOL" import "$WORK/view.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.lmdb.ssam" \
    --aln "$WORK/view.isam"
"$TOOL" sort --all "$WORK/view.isam"

seq_r1=$("$TOOL" view --no-header "$WORK/view.isam" \
    | awk '$1=="r1" { print $10; exit }')
[[ $seq_r1 == ACGTACGT ]]

[[ $("$TOOL" view -c "$WORK/view.isam") == 3 ]]

echo "lmdb store OK"
