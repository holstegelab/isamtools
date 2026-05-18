#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}
DATA_DIR=${DATA_DIR:-$ROOT/data}
CANON_BAM=${CANON_BAM:-$DATA_DIR/18R1995.merged.hifi.chm13.bam}
ALN_BAM=${ALN_BAM:-$DATA_DIR/18R1995.merged.hifi.hg38.bam}
N=${N:-10000}
SUBSET_N=${SUBSET_N:-500}
KEEP_WORK=${KEEP_WORK:-0}

[[ -x $TOOL ]] || { echo "build isamtools first (make)"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }
[[ -f $CANON_BAM ]] || { echo "missing canonical BAM: $CANON_BAM"; exit 1; }
[[ -f $ALN_BAM ]] || { echo "missing alignment BAM: $ALN_BAM"; exit 1; }

WORK=${WORK:-$(mktemp -d /tmp/isamtools-private-validate.XXXXXX)}
cleanup() {
    if [[ $KEEP_WORK != 0 ]]; then
        echo "kept workdir: $WORK" >&2
        return
    fi
    rm -rf "$WORK"
}
trap cleanup EXIT

pick_unique_qnames() {
    local bam=$1
    local n=$2
    local out=$3
    samtools view -F 2304 "$bam" \
        | awk -v limit="$n" '!seen[$1]++ && count < limit { print $1; ++count }' \
        > "$out"
    [[ -s $out ]] || { echo "no qnames selected from $bam"; exit 1; }
}

pick_unique_qnames "$CANON_BAM" "$N" "$WORK/qnames.txt"
samtools view -b -N "$WORK/qnames.txt" -o "$WORK/canonical.sub.bam" "$CANON_BAM"
samtools view -b -N "$WORK/qnames.txt" -o "$WORK/alignment.sub.bam" "$ALN_BAM"
samtools index "$WORK/canonical.sub.bam"
samtools index "$WORK/alignment.sub.bam"

"$TOOL" index "$WORK/canonical.sub.bam"
[[ -f "$WORK/canonical.sub.ssam" ]]
[[ ! -f "$WORK/canonical.sub.bam.ssam" ]]

"$TOOL" import "$WORK/alignment.sub.bam" \
    --store "$WORK/canonical.sub.bam" \
    --aln "$WORK/alignment.sub.isam" \
    --allow-missing
"$TOOL" sort --all "$WORK/alignment.sub.isam"

pick_unique_qnames "$WORK/alignment.sub.bam" "$SUBSET_N" "$WORK/hifi.txt"
"$TOOL" annotate "$WORK/canonical.sub.bam" --subset hifi --from "$WORK/hifi.txt"

subset_count=$($TOOL view -c --subset hifi "$WORK/alignment.sub.isam")
[[ $subset_count -gt 0 ]]

subset_region=$($TOOL view --no-header --subset hifi "$WORK/alignment.sub.isam" \
    | awk '$3 != "*" && !done { print $3 ":" $4 "-" $4; done = 1 }')
[[ -n $subset_region ]]
subset_region_count=$($TOOL view -c --subset hifi "$WORK/alignment.sub.isam" "$subset_region")
[[ $subset_region_count -gt 0 ]]

first_qname=$(head -n 1 "$WORK/hifi.txt")
canon_seq=$(samtools view "$WORK/canonical.sub.bam" \
    | awk -v q="$first_qname" '$1 == q && $2 < 256 && !done { print $10; done = 1 }')
view_seq=$($TOOL view --no-header --subset hifi "$WORK/alignment.sub.isam" \
    | awk -v q="$first_qname" '$1 == q && $2 < 256 && !done { print $10; done = 1 }')
[[ -n $canon_seq ]]
[[ -n $view_seq ]]
[[ $canon_seq == $view_seq ]]

cat > "$WORK/bad_store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qH	0	chrS	1	60	2H6M	*	0	0	ACGTAC	FFFFFF
SAM
samtools view -bS "$WORK/bad_store.sam" -o "$WORK/bad_store.bam"
if "$TOOL" index "$WORK/bad_store.bam" >/dev/null 2>"$WORK/bad_store.err"; then
    echo "expected hard-clipped canonical record to be rejected"
    exit 1
fi
grep -q "canonical store reads must keep full SEQ/QUAL" "$WORK/bad_store.err"

echo "private data validation OK"
echo "  workdir: $WORK"
echo "  canonical index: $WORK/canonical.sub.ssam"
echo "  alignment view:  $WORK/alignment.sub.isam"
echo "  subset count:    $subset_count"
echo "  subset region:   $subset_region ($subset_region_count records)"
