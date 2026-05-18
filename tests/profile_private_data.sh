#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}
DATA_DIR=${DATA_DIR:-$ROOT/data}
CANON_BAM=${CANON_BAM:-$DATA_DIR/18R1995.merged.hifi.chm13.bam}
ALN_BAM=${ALN_BAM:-$DATA_DIR/18R1995.merged.hifi.hg38.bam}
N=${N:-100000}
SUBSET_N=${SUBSET_N:-2000}
KEEP_WORK=${KEEP_WORK:-0}

[[ -x $TOOL ]] || { echo "build isamtools first (make)"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }
[[ -f $CANON_BAM ]] || { echo "missing canonical BAM: $CANON_BAM"; exit 1; }
[[ -f $ALN_BAM ]] || { echo "missing alignment BAM: $ALN_BAM"; exit 1; }

WORK=${WORK:-$(mktemp -d /tmp/isamtools-profile-private.XXXXXX)}
cleanup() {
    if [[ $KEEP_WORK != 0 ]]; then
        echo "kept workdir: $WORK" >&2
        return
    fi
    rm -rf "$WORK"
}
trap cleanup EXIT

run_timed() {
    local label=$1
    shift
    local start end ms
    start=$(date +%s%N)
    "$@"
    end=$(date +%s%N)
    ms=$(( (end - start) / 1000000 ))
    printf '%s\twall_ms=%d\n' "$label" "$ms"
}

pick_unique_qnames() {
    local bam=$1
    local n=$2
    local out=$3
    samtools view -F 2304 "$bam" \
        | awk -v limit="$n" '!seen[$1]++ && count < limit { print $1; ++count }' \
        > "$out"
    [[ -s $out ]] || { echo "no qnames selected from $bam"; exit 1; }
}

printf 'workdir\t%s\n' "$WORK"
printf 'reads\t%s\n' "$N"

run_timed sample_qnames \
    pick_unique_qnames "$CANON_BAM" "$N" "$WORK/qnames.txt"
run_timed subset_canonical \
    samtools view -b -N "$WORK/qnames.txt" -o "$WORK/canonical.sub.bam" "$CANON_BAM"
run_timed subset_alignment \
    samtools view -b -N "$WORK/qnames.txt" -o "$WORK/alignment.sub.bam" "$ALN_BAM"
run_timed samtools_index_canonical \
    samtools index "$WORK/canonical.sub.bam"
run_timed samtools_index_alignment \
    samtools index "$WORK/alignment.sub.bam"
run_timed isam_index \
    "$TOOL" index "$WORK/canonical.sub.bam"
run_timed isam_import \
    "$TOOL" import "$WORK/alignment.sub.bam" \
        --store "$WORK/canonical.sub.bam" \
        --aln "$WORK/alignment.sub.isam" \
        --allow-missing
run_timed isam_sort_all \
    "$TOOL" sort --all "$WORK/alignment.sub.isam"

pick_unique_qnames "$WORK/alignment.sub.bam" "$SUBSET_N" "$WORK/subset.txt"
run_timed annotate_store_subset \
    "$TOOL" annotate "$WORK/canonical.sub.bam" --subset bench --from "$WORK/subset.txt"
run_timed view_full \
    "$TOOL" view --no-header "$WORK/alignment.sub.isam" >/dev/null
run_timed view_subset \
    "$TOOL" view --no-header --subset bench "$WORK/alignment.sub.isam" >/dev/null

region=$($TOOL view --no-header "$WORK/alignment.sub.isam" \
    | awk '$3 != "*" && !done { print $3 ":" $4 "-" ($4 + 100000); done = 1 }')
if [[ -n $region ]]; then
    run_timed view_region \
        "$TOOL" view --no-header "$WORK/alignment.sub.isam" "$region" >/dev/null
    run_timed view_subset_region \
        "$TOOL" view --no-header --subset bench "$WORK/alignment.sub.isam" "$region" >/dev/null
    printf 'region\t%s\n' "$region"
fi

printf 'canonical_ssam\t%s\n' "$WORK/canonical.sub.ssam"
printf 'alignment_isam\t%s\n' "$WORK/alignment.sub.isam"
