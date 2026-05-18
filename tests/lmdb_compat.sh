#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)
TOOL=${TOOL:-$ROOT/build/isamtools}
BUILD_ROOT=$(cd "$(dirname "$TOOL")" && pwd)
HTSLIB_PREFIX=${HTSLIB_PREFIX:-$BUILD_ROOT/vendor/htslib}
LMDB_PREFIX=${LMDB_PREFIX:-$BUILD_ROOT/vendor/lmdb}
BZIP2_LIB=${BZIP2_LIB:-$BUILD_ROOT/vendor/bzip2/lib/libbz2.a}
HTSLIB_INC=$HTSLIB_PREFIX/include
HTSLIB_LIB=$HTSLIB_PREFIX/lib/libhts.a
LMDB_LIB=$LMDB_PREFIX/lib/liblmdb.a

if [[ -f $BZIP2_LIB ]]; then
    BZIP2_LINK=("$BZIP2_LIB")
else
    BZIP2_LINK=(-lbz2)
fi

[[ -x $TOOL ]] || { echo "built isamtools not found: $TOOL"; exit 1; }
command -v samtools >/dev/null || { echo "samtools is required"; exit 1; }
command -v cc >/dev/null || { echo "cc is required"; exit 1; }
[[ -f $HTSLIB_LIB ]] || { echo "local htslib static library not found: $HTSLIB_LIB"; exit 1; }
[[ -f $LMDB_LIB ]] || { echo "local LMDB static library not found: $LMDB_LIB"; exit 1; }

WORK=$(mktemp -d /tmp/isamtools-lmdb-compat.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

long_store_seq=$(printf '%*s' 500 '' | tr ' ' A)
long_store_qual=$(printf '%*s' 500 '' | tr ' ' F)
long_view_seq=$(printf '%*s' 500 '' | tr ' ' N)
long_view_qual=$(printf '%*s' 500 '' | tr ' ' '!')

{
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chrS\tLN:10000\n'
    printf '@RG\tID:rg1\tLB:libA\n'
    printf '@RG\tID:rg2\tLB:libB\n'
    printf 'qShort\t0\tchrS\t1\t60\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n'
    printf 'qLong\t0\tchrS\t10\t60\t500M\t*\t0\t0\t%s\t%s\n' \
        "$long_store_seq" "$long_store_qual"
    printf 'qOther\t0\tchrS\t800\t60\t8M\t*\t0\t0\tTTTTCCCC\tFFFFFFFF\n'
    printf 'qPair\t0\tchrS\t900\t60\t8M\t*\t0\t0\tAAAACCCC\tFFFFFFFF\n'
} > "$WORK/store.sam"

{
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chrV\tLN:10000\n'
    printf '@RG\tID:rg1\tLB:libA\n'
    printf '@RG\tID:rg2\tLB:libB\n'
    printf 'qShort\t0\tchrV\t100\t30\t8M\t*\t0\t0\tNNNNNNNN\t!!!!!!!!\tRG:Z:rg1\tHP:i:1\tZZ:Z:keep\n'
    printf 'qLong\t0\tchrV\t1\t40\t500M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\tHP:i:1\tZZ:Z:keep\n' \
        "$long_view_seq" "$long_view_qual"
    printf 'qOther\t0\tchrV\t800\t20\t8M\t*\t0\t0\tNNNNNNNN\t!!!!!!!!\tRG:Z:rg2\tHP:i:2\tZZ:Z:drop\n'
    printf 'qPair\t99\tchrV\t300\t50\t8M\t=\t900\t608\tNNNNNNNN\t!!!!!!!!\tRG:Z:rg1\tHP:i:3\tZZ:Z:pair\n'
    printf 'qPair\t147\tchrV\t900\t50\t8M\t=\t300\t-608\tNNNNNNNN\t!!!!!!!!\tRG:Z:rg1\tHP:i:3\tZZ:Z:pair\n'
} > "$WORK/aln.sam"

samtools view -bS "$WORK/store.sam" -o "$WORK/store.bam"

"$TOOL" index -o "$WORK/store.isam" "$WORK/store.bam"
test -f "$WORK/store.isam"
test -f "$WORK/store.isam-lock"

"$TOOL" import "$WORK/aln.sam" \
    --store "$WORK/store.bam" \
    --store-index "$WORK/store.isam" \
    --aln "$WORK/view.isam"
test -f "$WORK/view.isam"
test -f "$WORK/view.isam-lock"

"$TOOL" sort --interval "$WORK/view.isam"

printf 'qShort\nqLong\n' > "$WORK/subset.txt"
"$TOOL" annotate "$WORK/view.isam" --subset keep --from "$WORK/subset.txt"

[[ $("$TOOL" view -c "$WORK/view.isam" chrV:250-260) == 1 ]]
"$TOOL" view -c -o "$WORK/count.txt" "$WORK/view.isam" chrV:250-260
[[ $(cat "$WORK/count.txt") == 1 ]]

printf 'chrV\t249\t260\nchrV\t799\t807\n' > "$WORK/targets.bed"
[[ $("$TOOL" view -c -L "$WORK/targets.bed" "$WORK/view.isam") == 2 ]]
[[ $("$TOOL" view -c -L "$WORK/targets.bed" "$WORK/view.isam" chrV:250-260) == 1 ]]
[[ $("$TOOL" view -c "$WORK/view.isam" chrV:90-110 chrV:795-810) == 3 ]]
[[ $("$TOOL" view -c -L "$WORK/targets.bed" "$WORK/view.isam" chrV:1-1000 chrV:2000-3000) == 2 ]]
printf 'qShort\nqOther\nqMissing\n' > "$WORK/qnames.txt"
[[ $("$TOOL" view -c -N "$WORK/qnames.txt" "$WORK/view.isam") == 2 ]]
printf 'rg1\n' > "$WORK/rg.txt"
[[ $("$TOOL" view -c -r rg1 "$WORK/view.isam") == 4 ]]
[[ $("$TOOL" view -c -R "$WORK/rg.txt" "$WORK/view.isam") == 4 ]]
[[ $("$TOOL" view -c -l libB "$WORK/view.isam") == 1 ]]
[[ $("$TOOL" view -c -m 100 "$WORK/view.isam") == 1 ]]
[[ $("$TOOL" view -c -d ZZ:keep "$WORK/view.isam") == 2 ]]
printf 'pair\nkeep\n' > "$WORK/zz_values.txt"
[[ $("$TOOL" view -c -D ZZ:"$WORK/zz_values.txt" "$WORK/view.isam") == 4 ]]
[[ $("$TOOL" view -c -P "$WORK/view.isam" chrV:895-905) == 2 ]]
"$TOOL" view -H "$WORK/view.isam" | grep -q '^@SQ[[:space:]]SN:chrV'
if "$TOOL" view "$WORK/view.isam" | grep -q '^@'; then
    echo "default SAM view emitted header records"
    exit 1
fi
if "$TOOL" view --no-header "$WORK/view.isam" | grep -q '^@'; then
    echo "--no-header emitted header records"
    exit 1
fi
if ! "$TOOL" view -h "$WORK/view.isam" | grep -q '^@SQ[[:space:]]SN:chrV'; then
    echo "-h did not emit SAM header records"
    exit 1
fi
"$TOOL" view -O bam -o "$WORK/region.bam" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -c "$WORK/region.bam") == 1 ]]
"$TOOL" view -b -o "$WORK/region-b.bam" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -c "$WORK/region-b.bam") == 1 ]]
"$TOOL" view -1 -o "$WORK/region-fast.bam" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -c "$WORK/region-fast.bam") == 1 ]]
"$TOOL" view -u -o "$WORK/region-uncompressed.bam" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -c "$WORK/region-uncompressed.bam") == 1 ]]
"$TOOL" view -O bam,level=1 -o "$WORK/region-level.bam" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -c "$WORK/region-level.bam") == 1 ]]

{
    printf '>chrV\n'
    printf '%*s\n' 10000 '' | tr ' ' A
} > "$WORK/ref.fa"
samtools faidx "$WORK/ref.fa"
"$TOOL" view -C -T "$WORK/ref.fa" -o "$WORK/region.cram" "$WORK/view.isam" chrV:250-260
[[ $(samtools view -T "$WORK/ref.fa" -c "$WORK/region.cram") == 1 ]]

cat > "$WORK/c_api_lmdb.c" <<'C'
#include <stdio.h>
#include <string.h>
#include <htslib/sam.h>
#include <isam/isam.h>

static int fail(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    return 1;
}

static int count_itr(isamFile *fp, isamItr *itr) {
    bam1_t *b = bam_init1();
    int n = 0;
    int rc;
    while ((rc = isam_itr_next(fp, itr, b)) >= 0) ++n;
    bam_destroy1(b);
    return rc == ISAM_EOF ? n : rc;
}

int main(int argc, char **argv) {
    if (argc != 2) return fail("usage: c_api_lmdb <view.isam>");

    isamFile *fp = isam_open(argv[1], "r");
    if (!fp) return fail("isam_open failed");
    sam_hdr_t *hdr = isam_hdr_read(fp);
    if (!hdr) return fail("isam_hdr_read failed");
    bam1_t *b = bam_init1();
    int seq_count = 0;
    int qlong_count = 0;
    int rc;
    while ((rc = isam_read1(fp, hdr, b)) >= 0) {
        ++seq_count;
        if (strcmp(bam_get_qname(b), "qLong") == 0) ++qlong_count;
    }
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    isam_close(fp);
    if (rc != ISAM_EOF) return fail("isam_read1 returned an error");
    if (seq_count != 5 || qlong_count != 1) return fail("unexpected sequential counts");

    fp = isam_open(argv[1], "r");
    if (!fp) return fail("second isam_open failed");
    hdr = isam_hdr_read(fp);
    isamIndex *idx = isam_index_load(argv[1]);
    if (!idx) return fail("isam_index_load failed");

    isamItr *region = isam_itr_querys(idx, hdr, "chrV:250-260");
    if (!region) return fail("isam_itr_querys failed");
    int region_count = count_itr(fp, region);
    isam_itr_destroy(region);
    if (region_count != 1) return fail("unexpected region count");

    isamItr *name = isam_itr_query_name(idx, "qLong");
    if (!name) return fail("isam_itr_query_name failed");
    int name_count = count_itr(fp, name);
    isam_itr_destroy(name);
    if (name_count != 1) return fail("unexpected qname count");

    isamItr *subset = isam_itr_query_subset(idx, "keep");
    if (!subset) return fail("isam_itr_query_subset failed");
    int subset_count = count_itr(fp, subset);
    isam_itr_destroy(subset);
    if (subset_count != 2) return fail("unexpected subset count");

    isamItr *region_subset = isam_itr_query_region_subset(idx, hdr,
                                                          "chrV:250-260",
                                                          "keep");
    if (!region_subset) return fail("isam_itr_query_region_subset failed");
    int region_subset_count = count_itr(fp, region_subset);
    isam_itr_destroy(region_subset);
    if (region_subset_count != 1) return fail("unexpected region subset count");

    isam_index_destroy(idx);
    sam_hdr_destroy(hdr);
    isam_close(fp);
    puts("c api lmdb OK");
    return 0;
}
C

cc -no-pie -I"$ROOT/include" -I"$HTSLIB_INC" \
    "$WORK/c_api_lmdb.c" "$BUILD_ROOT/libisam.a" \
    "$LMDB_LIB" "$HTSLIB_LIB" "${BZIP2_LINK[@]}" -lz -llzma -lm -lpthread -lzstd -o "$WORK/c_api_lmdb"
"$WORK/c_api_lmdb" "$WORK/view.isam"

echo "lmdb compat OK"
