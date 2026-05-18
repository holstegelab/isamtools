#!/usr/bin/env bash
# Smoke test for the full isamtools pipeline using a small chrM subset.
# Builds a canonical chm13 store + an hg38 alignment-view, sorts, and
# asserts that reconstructed SEQ for one read matches the canonical SEQ.
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}
TD=$ROOT/test_data
DATA=$ROOT/data

[[ -x $TOOL ]] || { echo "build isamtools first (make)"; exit 1; }

mkdir -p "$TD"

# 1. Build small subsets if they don't exist.
if [[ ! -s $TD/chm13.sub.bam ]]; then
    samtools view -b -o "$TD/chm13.sub.bam" "$DATA/18R1995.merged.hifi.chm13.bam" chrM
    samtools index "$TD/chm13.sub.bam"
fi
samtools view -F 2304 "$TD/chm13.sub.bam" | cut -f1 | sort -u > "$TD/qnames.txt"
if [[ ! -s $TD/hg38.sub.bam ]]; then
    samtools view -b -N "$TD/qnames.txt" -o "$TD/hg38.sub.bam" \
        "$DATA/18R1995.merged.hifi.hg38.bam" chrM
fi

# 2. Pipeline.
rm -f "$TD"/*.isam "$TD"/*.ssam
"$TOOL" index  "$TD/chm13.sub.bam"
"$TOOL" import "$TD/hg38.sub.bam" \
    --store "$TD/chm13.sub.bam" \
    --aln   "$TD/hg38.sub.isam" \
    --allow-missing
"$TOOL" sort   --all "$TD/hg38.sub.isam"

# 3. Sanity: reconstructed SEQ md5 must equal canonical SEQ md5.
qname=$(head -1 "$TD/qnames.txt")
recon=$("$TOOL" view "$TD/hg38.sub.isam" \
        | awk -v q="$qname" '!seen && $1==q && $2 < 256 {print $10; seen=1}' | md5sum)
canon=$(samtools view "$TD/chm13.sub.bam" \
        | awk -v q="$qname" '!seen && $1==q && $2 < 256 {print $10; seen=1}' | md5sum)

if [[ "$recon" != "$canon" ]]; then
    echo "MISMATCH: reconstructed seq differs from canonical"
    echo "  recon: $recon"
    echo "  canon: $canon"
    exit 1
fi
echo "smoke OK: reconstructed SEQ matches canonical ($recon)"
