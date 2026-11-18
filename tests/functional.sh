#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}

[[ -x $TOOL ]] || { echo "build isamtools first (make)"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-functional.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/store_a.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF	RG:Z:canon
qB	0	chrS	20	50	8M	*	0	0	TTTTCCCC	FFFFFFFF	RG:Z:canon
qC	4	*	0	0	*	*	0	0	GGGGTTTT	FFFFFFFF	RG:Z:canon
qD	0	chrS	40	50	8M	*	0	0	ACGTACGT	FFFFFFFF	RG:Z:canon
SAM

cat > "$WORK/store_b.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qB	0	chrS	20	50	8M	*	0	0	TTTTCCCC	FFFFFFFF	RG:Z:canon
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF	RG:Z:canon
qC	4	*	0	0	*	*	0	0	GGGGTTTT	FFFFFFFF	RG:Z:canon
qD	0	chrS	40	50	8M	*	0	0	ACGTACGT	FFFFFFFF	RG:Z:canon
SAM

cat > "$WORK/aln.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qA	0	chrV	101	30	8M	*	0	0	NNNNNNNN	!!!!!!!!	NM:i:0
qB	16	chrV	200	20	8M	*	0	0	NNNNNNNN	!!!!!!!!	NM:i:1
qA	2048	chrV	150	10	4M4S	*	0	0	NNNNNNNN	!!!!!!!!	SA:Z:chrV,101,+,8M,30,0;
qC	4	*	0	0	*	*	0	0	NNNNNNNN	!!!!!!!!
qD	0	chrV	300	30	2H6M	*	0	0	NNNNNN	!!!!!!
SAM

cat > "$WORK/missing.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qZ	0	chrV	1	60	8M	*	0	0	AAAAAAAA	FFFFFFFF
SAM

cat > "$WORK/fallback_store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qP	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
qSupp	2048	chrS	10	60	8M	*	0	0	TTTTCCCC	FFFFFFFF	SA:Z:chrS,10,+,8M,60,0;
qHard	2048	chrS	20	60	2H6M	*	0	0	GGGGTT	FFFFFF	SA:Z:chrS,20,+,2H6M,60,0;
SAM

cat > "$WORK/fallback_import.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qSupp	0	chrV	500	60	8M	*	0	0	NNNNNNNN	!!!!!!!!
SAM

cat > "$WORK/fallback_missing.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrV	LN:1000
qHard	0	chrV	600	60	6M	*	0	0	NNNNNN	!!!!!!
SAM

samtools view -bS "$WORK/store_a.sam" -o "$WORK/store_a.bam"
samtools view -bS "$WORK/store_b.sam" -o "$WORK/store_b.bam"
samtools view -bS "$WORK/fallback_store.sam" -o "$WORK/fallback_store.bam"

"$TOOL" index -@ 2 "$WORK/store_a.bam" 2>"$WORK/store_a.index.err"
grep -q "input decompression threads: 2" "$WORK/store_a.index.err"
[[ -f "$WORK/store_a.ssam" ]]
[[ ! -f "$WORK/store_a.bam.ssam" ]]

"$TOOL" check "$WORK/store_a.bam" --store "$WORK/store_a.bam" 2>"$WORK/check_store.err"
grep -q "4 primary records match store" "$WORK/check_store.err"

awk 'BEGIN { print ">chrS"; for (i = 0; i < 1000; ++i) printf "A"; print "" }' \
    > "$WORK/reference.fa"
samtools faidx "$WORK/reference.fa"
samtools sort -@ 2 -O CRAM --reference "$WORK/reference.fa" \
    -o "$WORK/store_a.cram" "$WORK/store_a.bam"
samtools index "$WORK/store_a.cram"
"$TOOL" index -@ 2 -T "$WORK/reference.fa" \
    -o "$WORK/store_a.cram.ssam" "$WORK/store_a.cram" \
    2>"$WORK/store_a.cram.index.err"
grep -q "input decompression threads: 2" "$WORK/store_a.cram.index.err"
"$TOOL" import -@ 2 --fast-verify "$WORK/store_a.bam" \
    --store "$WORK/store_a.cram" \
    --store-index "$WORK/store_a.cram.ssam" \
    --aln "$WORK/store_a.cram.isam" 2>"$WORK/store_a.cram.import.err"
"$TOOL" check "$WORK/store_a.bam" --view "$WORK/store_a.cram.isam" \
    2>"$WORK/store_a.cram.check.err"
grep -q "records match reconstructed view" "$WORK/store_a.cram.check.err"

# A view may revisit CRAM containers in an interleaved order.  Check that a
# bounded batch groups those payload fetches by container while retaining the
# original record order for reconstruction and comparison.
awk 'BEGIN {
    print ">chrB"
    for (i = 0; i < 2000000; ++i) printf "A"
    print ""
}' > "$WORK/batch_reference.fa"
samtools faidx "$WORK/batch_reference.fa"
awk 'BEGIN {
    OFS = "\t"
    print "@HD", "VN:1.6", "SO:coordinate"
    print "@SQ", "SN:chrB", "LN:2000000"
    for (i = 0; i < 10000; ++i) {
        seq = seq "A"
        qual = qual "F"
    }
    for (container = 0; container < 80; ++container)
        for (within = 0; within < 2; ++within) {
            qname = sprintf("qBatch%03d_%d", container, within)
            pos = container + 1
            print qname, 0, "chrB", pos, 60, "10000M", "*", 0, 0, seq, qual
        }
}' > "$WORK/batch_store.sam"
awk 'BEGIN {
    OFS = "\t"
    print "@HD", "VN:1.6", "SO:unsorted"
    print "@SQ", "SN:chrB", "LN:2000000"
    for (i = 0; i < 10000; ++i) {
        seq = seq "A"
        qual = qual "F"
    }
    for (within = 0; within < 2; ++within)
        for (container = 0; container < 80; ++container) {
            qname = sprintf("qBatch%03d_%d", container, within)
            pos = container + 1
            print qname, 0, "chrB", pos, 60, "10000M", "*", 0, 0, seq, qual
        }
}' > "$WORK/batch_input.sam"
samtools view -bS "$WORK/batch_input.sam" -o "$WORK/batch_input.bam"
samtools view -C -T "$WORK/batch_reference.fa" \
    --output-fmt-option seqs_per_slice=32 \
    --output-fmt-option slices_per_container=1 \
    -o "$WORK/batch_store.cram" "$WORK/batch_store.sam"
samtools index "$WORK/batch_store.cram"
"$TOOL" index -T "$WORK/batch_reference.fa" \
    -o "$WORK/batch_store.ssam" "$WORK/batch_store.cram" \
    2>"$WORK/batch_store.index.err"
env CRAM_REFERENCE="$WORK/batch_reference.fa" "$TOOL" import --fast-verify \
    "$WORK/batch_input.bam" \
    --store "$WORK/batch_store.cram" \
    --store-index "$WORK/batch_store.ssam" \
    --aln "$WORK/batch_view.isam" 2>"$WORK/batch_import.err"
env CRAM_REFERENCE="$WORK/batch_reference.fa" ISAM_FETCH_STATS=1 \
    ISAM_SLICE_CACHE_MB=1 "$TOOL" check --batch-records 1 \
    "$WORK/batch_input.bam" --view "$WORK/batch_view.isam" \
    --store "$WORK/batch_store.cram" \
    --store-index "$WORK/batch_store.ssam" 2>"$WORK/batch_single.err"
env CRAM_REFERENCE="$WORK/batch_reference.fa" ISAM_FETCH_STATS=1 \
    ISAM_SLICE_CACHE_MB=1 "$TOOL" check --batch-records 160 \
    "$WORK/batch_input.bam" --view "$WORK/batch_view.isam" \
    --store "$WORK/batch_store.cram" \
    --store-index "$WORK/batch_store.ssam" 2>"$WORK/batch_grouped.err"
grep -q "160 records match reconstructed view" "$WORK/batch_single.err"
grep -q "160 records match reconstructed view" "$WORK/batch_grouped.err"
single_decodes=$(sed -n 's/.*miss\/decode=\([0-9][0-9]*\).*/\1/p' \
    "$WORK/batch_single.err" | tail -1)
grouped_decodes=$(sed -n 's/.*miss\/decode=\([0-9][0-9]*\).*/\1/p' \
    "$WORK/batch_grouped.err" | tail -1)
[[ $single_decodes =~ ^[0-9]+$ && $grouped_decodes =~ ^[0-9]+$ ]]
(( grouped_decodes < single_decodes )) || {
    echo "expected batched check to reduce CRAM container decodes ($single_decodes vs $grouped_decodes)"
    exit 1
}
echo "batched check reduced CRAM container decodes: $single_decodes -> $grouped_decodes"

# Output compression threads must not change reconstruction or coordinate
# ordering when CRAM payloads revisit several container locators.
"$TOOL" sort --all "$WORK/batch_view.isam"
for workers in 1 4; do
    env ISAM_SLICE_CACHE_MB=1 "$TOOL" view -@ "$workers" -b \
        -T "$WORK/batch_reference.fa" --store "$WORK/batch_store.cram" \
        --store-index "$WORK/batch_store.ssam" \
        -o "$WORK/batch_view.t${workers}.bam" "$WORK/batch_view.isam"
    samtools quickcheck "$WORK/batch_view.t${workers}.bam"
    samtools index "$WORK/batch_view.t${workers}.bam"
    samtools view "$WORK/batch_view.t${workers}.bam" | LC_ALL=C sort \
        > "$WORK/batch_view.t${workers}.sam"
done
cmp "$WORK/batch_view.t1.sam" "$WORK/batch_view.t4.sam"
[[ $(samtools view -c "$WORK/batch_view.t4.bam") == 160 ]]

# CRAM is allowed to decode explicit =/X CIGAR operations as M, and HTSlib
# may choose different BAM integer widths for the same SAM i value.  Neither
# representation change should make a logically identical record fail check.
cat > "$WORK/cram_normalized.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qNorm	0	chrS	101	60	4=1X3=	*	0	0	AAAATAAA	FFFFFFFF	NM:i:1	A1:i:-1	A2:i:255	A3:i:65535	A4:i:2147483647
qUnmap255	4	*	0	255	*	*	0	0	ACGT	FFFF
SAM
samtools view -bS "$WORK/cram_normalized.sam" -o "$WORK/cram_normalized.bam"
samtools view -C -T "$WORK/reference.fa" \
    -o "$WORK/cram_normalized.cram" "$WORK/cram_normalized.bam"
samtools index "$WORK/cram_normalized.cram"
"$TOOL" index -T "$WORK/reference.fa" "$WORK/cram_normalized.cram" \
    2>"$WORK/cram_normalized.index.err"
"$TOOL" check "$WORK/cram_normalized.bam" \
    --store "$WORK/cram_normalized.cram" \
    2>"$WORK/cram_normalized.store_check.err"
grep -q "2 primary records match store" "$WORK/cram_normalized.store_check.err"
"$TOOL" import --fast-verify "$WORK/cram_normalized.bam" \
    --store "$WORK/cram_normalized.cram" \
    --aln "$WORK/cram_normalized.isam" \
    2>"$WORK/cram_normalized.import.err"
"$TOOL" check "$WORK/cram_normalized.bam" \
    --view "$WORK/cram_normalized.isam" \
    2>"$WORK/cram_normalized.view_check.err"
grep -q "2 records match reconstructed view" "$WORK/cram_normalized.view_check.err"

cp "$WORK/cram_normalized.bam" "$WORK/semantic_expected.bam"
sed 's/4=1X3=/3=2X3=/' "$WORK/cram_normalized.sam" \
    > "$WORK/cigar_mismatch.sam"
samtools view -bS "$WORK/cigar_mismatch.sam" -o "$WORK/cigar_mismatch.bam"
samtools index "$WORK/cigar_mismatch.bam"
"$TOOL" index "$WORK/cigar_mismatch.bam" \
    2>"$WORK/cigar_mismatch.index.err"
if "$TOOL" check "$WORK/semantic_expected.bam" \
        --store "$WORK/cigar_mismatch.bam" \
        2>"$WORK/cigar_mismatch.check.err"; then
    echo "expected check to reject different explicit =/X CIGAR operations"
    exit 1
fi
grep -q "differs from store" "$WORK/cigar_mismatch.check.err"

sed 's/A2:i:255/A2:i:254/' "$WORK/cram_normalized.sam" \
    > "$WORK/integer_mismatch.sam"
samtools view -bS "$WORK/integer_mismatch.sam" -o "$WORK/integer_mismatch.bam"
samtools index "$WORK/integer_mismatch.bam"
"$TOOL" index "$WORK/integer_mismatch.bam" \
    2>"$WORK/integer_mismatch.index.err"
if "$TOOL" check "$WORK/semantic_expected.bam" \
        --store "$WORK/integer_mismatch.bam" \
        2>"$WORK/integer_mismatch.check.err"; then
    echo "expected check to reject different integer aux values"
    exit 1
fi
grep -q "differs from store" "$WORK/integer_mismatch.check.err"

sed 's/	60	4=1X3=/	59	4=1X3=/' "$WORK/cram_normalized.sam" \
    > "$WORK/mapq_mismatch.sam"
samtools view -bS "$WORK/mapq_mismatch.sam" -o "$WORK/mapq_mismatch.bam"
samtools index "$WORK/mapq_mismatch.bam"
"$TOOL" index "$WORK/mapq_mismatch.bam" \
    2>"$WORK/mapq_mismatch.index.err"
if "$TOOL" check "$WORK/semantic_expected.bam" \
        --store "$WORK/mapq_mismatch.bam" \
        2>"$WORK/mapq_mismatch.check.err"; then
    echo "expected check to reject different mapped MAPQ values"
    exit 1
fi
grep -q "differs from store" "$WORK/mapq_mismatch.check.err"

if "$TOOL" check "$WORK/aln.sam" --store "$WORK/store_a.bam" >/dev/null 2>"$WORK/check_mismatch.err"; then
    echo "expected import --check to fail when input primary records differ from store"
    exit 1
fi
grep -q "differs from store" "$WORK/check_mismatch.err"

"$TOOL" import --check --threads 2 "$WORK/store_a.bam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/check_import.isam" 2>"$WORK/check_import.err"
grep -q "worker threads: 2 (1 input decompression, 1 record compression)" "$WORK/check_import.err"
grep -q "records match reconstructed view" "$WORK/check_import.err"
"$TOOL" check "$WORK/store_a.bam" --view "$WORK/check_import.isam" 2>"$WORK/check_view.err"
grep -q "records match reconstructed view" "$WORK/check_view.err"

cat "$WORK/store_a.bam" | env ISAM_FETCH_STATS=1 "$TOOL" import -@ 2 --fast-verify - \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/threaded_stdin.isam" 2>"$WORK/threaded_stdin.err"
grep -q "worker threads: 2 (1 input decompression, 1 record compression)" "$WORK/threaded_stdin.err"
grep -q "4 fingerprint-verified" "$WORK/threaded_stdin.err"
grep -q "voffset=0" "$WORK/threaded_stdin.err"
[[ $($TOOL view -c "$WORK/threaded_stdin.isam") == 4 ]]

# Exercise the ordered compression queue after its 1000-record dictionary
# training prefix. Repeating the canonical records keeps this fixture compact
# while still checking threaded import, queue flush, order, and reconstruction.
samtools view -h "$WORK/store_a.bam" | awk '
    /^@/ { print; next }
    { records[++n] = $0 }
    END {
        for (round = 0; round < 275; ++round)
            for (i = 1; i <= n; ++i)
                print records[i]
    }
' | samtools view -bS - -o "$WORK/threaded_large.bam"
"$TOOL" import -@ 2 --check --fast-verify "$WORK/threaded_large.bam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/threaded_large.isam" 2>"$WORK/threaded_large.err"
grep -q "worker threads: 2 (1 input decompression, 1 record compression)" \
    "$WORK/threaded_large.err"
grep -q "trained aln_by_id zstd dictionary from 1000 records" \
    "$WORK/threaded_large.err"
grep -q "1100 records match reconstructed view" "$WORK/threaded_large.err"
[[ $($TOOL view -c "$WORK/threaded_large.isam") == 1100 ]]

cat > "$WORK/fingerprint_mismatch.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qA	0	chrS	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF	RG:Z:alternate
SAM
env ISAM_FETCH_STATS=1 "$TOOL" import --fast-verify \
    "$WORK/fingerprint_mismatch.sam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/fingerprint_mismatch.isam" 2>"$WORK/fingerprint_mismatch.err"
grep -q "1 fingerprint-verified (0 fields matched, 1 differed/absent)" \
    "$WORK/fingerprint_mismatch.err"
grep -q "voffset=0" "$WORK/fingerprint_mismatch.err"
"$TOOL" view "$WORK/fingerprint_mismatch.isam" | grep -q $'\tRG:Z:alternate'

"$TOOL" index "$WORK/fallback_store.bam" >/dev/null 2>"$WORK/fallback_store.err"
grep -q "used supplementary alignments as canonical records for 1 reads without a primary alignment" "$WORK/fallback_store.err"
grep -q "1 reads had only hard-clipped supplementary alignments and were not indexed" "$WORK/fallback_store.err"

"$TOOL" import "$WORK/fallback_import.sam" \
    --store "$WORK/fallback_store.bam" \
    --aln "$WORK/fallback_view.isam"
[[ $($TOOL view -c "$WORK/fallback_view.isam") == 1 ]]

if "$TOOL" import "$WORK/fallback_missing.sam" \
        --store "$WORK/fallback_store.bam" \
        --aln "$WORK/fallback_missing.isam" >/dev/null 2>"$WORK/fallback_missing.err"; then
    echo "expected import to fail when only hard-clipped supplementary canonical record exists"
    exit 1
fi
grep -q "falls back to non-hard-clipped supplementary alignments" "$WORK/fallback_missing.err"

cat > "$WORK/bad_store.sam" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chrS	LN:1000
qH	0	chrS	1	60	2H6M	*	0	0	ACGTAC	FFFFFF
SAM
samtools view -bS "$WORK/bad_store.sam" -o "$WORK/bad_store.bam"

if "$TOOL" index "$WORK/bad_store.bam" >/dev/null 2>"$WORK/bad_store.err"; then
    echo "expected index to fail on hard-clipped canonical primary read"
    exit 1
fi
grep -q "canonical store reads must keep full SEQ/QUAL" "$WORK/bad_store.err"

if "$TOOL" import "$WORK/missing.sam" \
        --store "$WORK/store_a.bam" \
        --aln "$WORK/missing.isam" >/dev/null 2>"$WORK/missing.err"; then
    echo "expected import to fail on missing qname"
    exit 1
fi
grep -q "falls back to non-hard-clipped supplementary alignments" "$WORK/missing.err"

[[ -f "$WORK/missing.isam" ]]
"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/missing.isam"
[[ $($TOOL view -c "$WORK/missing.isam") == 5 ]]

cp "$WORK/aln.sam" "$WORK/default_out.sam"
"$TOOL" import "$WORK/default_out.sam" --store "$WORK/store_a.bam"
[[ -f "$WORK/default_out.isam" ]]
[[ $($TOOL view -c "$WORK/default_out.isam") == 5 ]]

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/view.isam"

[[ $("$TOOL" view --no-header "$WORK/view.isam" | awk 'END { print NR }') == 5 ]]
[[ $("$TOOL" view --no-header --order name "$WORK/view.isam" | awk 'END { print NR }') == 5 ]]

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/name_only.isam"
"$TOOL" sort -n "$WORK/name_only.isam"
[[ $("$TOOL" view --no-header --order name "$WORK/name_only.isam" | awk 'END { print NR }') == 5 ]]

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/interval_only.isam"
"$TOOL" sort --interval "$WORK/interval_only.isam"
[[ $("$TOOL" view --no-header "$WORK/interval_only.isam" chrV:145-155 | awk 'END { print NR }') == 1 ]]

"$TOOL" sort --all "$WORK/view.isam"

[[ $("$TOOL" view -c "$WORK/view.isam") == 5 ]]
[[ $("$TOOL" view --no-header --order name "$WORK/view.isam" | awk 'END { print NR }') == 5 ]]

seq_a=$("$TOOL" view --no-header "$WORK/view.isam" \
    | awk '$1=="qA" && $2==0 { print $10; exit }')
[[ $seq_a == ACGTACGT ]]

region_count=$("$TOOL" view --no-header "$WORK/view.isam" chrV:145-155 \
    | awk 'END { print NR }')
[[ $region_count == 1 ]]

hard_clip_seq=$("$TOOL" view --no-header "$WORK/view.isam" \
    | awk '$1=="qD" { print $10; exit }')
[[ $hard_clip_seq == GTACGT ]]
[[ $("$TOOL" view -c -f 2048 "$WORK/view.isam") == 1 ]]
[[ $("$TOOL" view -c -F 2048 "$WORK/view.isam") == 4 ]]

printf 'qA\tHP:i:1\n' > "$WORK/tags.tsv"
printf 'qA\n' > "$WORK/subset.txt"
"$TOOL" annotate "$WORK/view.isam" --from "$WORK/tags.tsv"
"$TOOL" annotate "$WORK/view.isam" --subset keep --from "$WORK/subset.txt"
[[ $("$TOOL" view -c --subset keep "$WORK/view.isam") == 2 ]]

printf 'qA\tST:Z:store-a\n' > "$WORK/store_tags.tsv"
"$TOOL" annotate "$WORK/store_a.bam" --from "$WORK/store_tags.tsv"
"$TOOL" annotate "$WORK/store_a.bam" --subset keep_store --from "$WORK/subset.txt"
"$TOOL" annotate "$WORK/store_a.bam" --subset all_from_bam --from-bam "$WORK/store_a.bam"

[[ $($TOOL view -c --subset keep_store "$WORK/view.isam") == 2 ]]
[[ $($TOOL view -c --subset all_from_bam "$WORK/view.isam") == 5 ]]

line_with_store_tag=$($TOOL view --no-header "$WORK/view.isam" | awk '$1=="qA" && $2==0 { print; exit }')
grep -q $'ST:Z:store-a' <<<"$line_with_store_tag"

"$TOOL" index -o "$WORK/store_b.ssam" "$WORK/store_b.bam"

"$TOOL" rehint "$WORK/view.isam" \
    --store "$WORK/store_b.bam" \
    --store-index "$WORK/store_b.ssam"

[[ $("$TOOL" view -c --subset keep "$WORK/view.isam") == 2 ]]
[[ $("$TOOL" view -c --subset keep "$WORK/view.isam" chrV:145-155) == 1 ]]

seq_a_rehint=$("$TOOL" view --no-header "$WORK/view.isam" \
    | awk '$1=="qA" && $2==0 { print $10; exit }')
[[ $seq_a_rehint == ACGTACGT ]]

echo "functional OK"
