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

"$TOOL" index "$WORK/store_a.bam"
[[ -f "$WORK/store_a.ssam" ]]
[[ ! -f "$WORK/store_a.bam.ssam" ]]

"$TOOL" check "$WORK/store_a.bam" --store "$WORK/store_a.bam" 2>"$WORK/check_store.err"
grep -q "4 primary records match store" "$WORK/check_store.err"

if "$TOOL" check "$WORK/aln.sam" --store "$WORK/store_a.bam" >/dev/null 2>"$WORK/check_mismatch.err"; then
    echo "expected import --check to fail when input primary records differ from store"
    exit 1
fi
grep -q "differs from store" "$WORK/check_mismatch.err"

"$TOOL" import --check "$WORK/store_a.bam" \
    --store "$WORK/store_a.bam" \
    --aln "$WORK/check_import.isam" 2>"$WORK/check_import.err"
grep -q "records match reconstructed view" "$WORK/check_import.err"
"$TOOL" check "$WORK/store_a.bam" --view "$WORK/check_import.isam" 2>"$WORK/check_view.err"
grep -q "records match reconstructed view" "$WORK/check_view.err"

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
