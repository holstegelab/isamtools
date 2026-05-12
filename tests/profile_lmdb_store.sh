#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}
N=${N:-5000}

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-profile-lmdb.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

run_timed() {
    local label=$1
    shift
    local start end ms
    start=$(date +%s%N)
    "$@"
    end=$(date +%s%N)
    ms=$(( (end - start) / 1000000 ))
    printf '%s wall_ms=%d\n' "$label" "$ms" >&2
}

awk -v n="$N" 'BEGIN {
    print "@HD\tVN:1.6\tSO:unsorted";
    print "@SQ\tSN:chrS\tLN:" (n * 5 + 1000);
    seq = ""; qual = "";
    for (j = 0; j < 100; ++j) { seq = seq "A"; qual = qual "F"; }
    for (i = 1; i <= n; ++i)
        printf "q%06d\t0\tchrS\t%d\t60\t100M\t*\t0\t0\t%s\t%s\n",
               i, i * 3, seq, qual;
}' > "$WORK/store.sam"

awk -v n="$N" 'BEGIN {
    print "@HD\tVN:1.6\tSO:unsorted";
    print "@SQ\tSN:chrV\tLN:" (n * 5 + 1000);
    seq = ""; qual = "";
    for (j = 0; j < 100; ++j) { seq = seq "N"; qual = qual "!"; }
    for (i = 1; i <= n; ++i)
        printf "q%06d\t0\tchrV\t%d\t40\t100M\t*\t0\t0\t%s\t%s\tNM:i:0\n",
               i, i * 3 + 1, seq, qual;
}' > "$WORK/aln.sam"

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

run_timed index_lmdb \
    "$TOOL" index -o "$WORK/store.lmdb.isam" "$WORK/store.bam"
run_timed import_lmdb_store \
    "$TOOL" import "$WORK/aln.sam" \
        --store "$WORK/store.bam" \
        --store-index "$WORK/store.lmdb.isam" \
        --aln "$WORK/view.isam"
run_timed sort_all \
    "$TOOL" sort --all "$WORK/view.isam"
run_timed view_full \
    "$TOOL" view --no-header "$WORK/view.isam" >/dev/null
run_timed view_region \
    "$TOOL" view --no-header "$WORK/view.isam" chrV:7500-8500 >/dev/null
