#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=$(pwd)

TOOL=${TOOL:-$ROOT/build/isamtools}
SAMTOOLS=${SAMTOOLS:-samtools}
CC=${CC:-cc}
SYNTHETIC=${SYNTHETIC:-1}
PREPARE=${PREPARE:-1}
KEEP_WORK=${KEEP_WORK:-0}
REPEATS=${REPEATS:-3}
N=${N:-20000}
SAMPLE_N=${SAMPLE_N:-0}
QNAME_N=${QNAME_N:-10000}
REGION_COUNT=${REGION_COUNT:-32}
REGION_WINDOW=${REGION_WINDOW:-10000}
TARGET_COUNT=${TARGET_COUNT:-128}
SPAN_DIST=${SPAN_DIST:-10000}
THREADS=${THREADS:-0}
RUN_C_API=${RUN_C_API:-1}
RUN_PERF=${RUN_PERF:-0}
ALLOW_MISSING=${ALLOW_MISSING:-0}
TIME_BIN=${TIME_BIN:-}
DATA_DIR=${DATA_DIR:-$ROOT/data}
CANON_BAM=${CANON_BAM:-$DATA_DIR/18R1995.merged.hifi.chm13.bam}
ALN_BAM=${ALN_BAM:-$DATA_DIR/18R1995.merged.hifi.hg38.bam}
STORE_SSAM=${STORE_SSAM:-}
ALN_ISAM=${ALN_ISAM:-}

[[ -x $TOOL ]] || { echo "build isamtools first (make)" >&2; exit 1; }
command -v "$SAMTOOLS" >/dev/null || { echo "samtools is required" >&2; exit 1; }
command -v "$CC" >/dev/null || { echo "cc is required" >&2; exit 1; }
[[ $REPEATS =~ ^[0-9]+$ && $REPEATS -gt 0 ]] || { echo "REPEATS must be a positive integer" >&2; exit 1; }
if [[ -z $TIME_BIN ]]; then
    if [[ -x /usr/bin/time ]]; then
        TIME_BIN=/usr/bin/time
    elif command -v gtime >/dev/null; then
        TIME_BIN=$(command -v gtime)
    fi
fi

WORK=${WORK:-$(mktemp -d /tmp/isamtools-profile-comprehensive.XXXXXX)}
REPORT_DIR=${REPORT_DIR:-$WORK/report}
mkdir -p "$WORK" "$REPORT_DIR"

cleanup() {
    if [[ $KEEP_WORK != 0 ]]; then
        echo "kept workdir: $WORK" >&2
        return
    fi
    rm -rf "$WORK"
}
trap cleanup EXIT

RAW_TSV=$REPORT_DIR/raw.tsv
SUMMARY_TSV=$REPORT_DIR/summary.tsv
REPORT_MD=$REPORT_DIR/report.md
COMMANDS_TSV=$REPORT_DIR/commands.tsv
ENV_TXT=$REPORT_DIR/environment.txt

printf 'case_id\tgroup\tengine\tscenario\trepeat\telapsed_sec\tuser_sec\tsys_sec\tmax_rss_kb\tstdout_bytes\tstderr_bytes\texit_status\tcommand\n' > "$RAW_TSV"
printf 'case_id\tgroup\tengine\tscenario\tcommand\n' > "$COMMANDS_TSV"

log() {
    printf '[profile] %s\n' "$*" >&2
}

shell_quote() {
    printf "%q" "$1"
}

join_quoted_file() {
    local file=$1
    local out=()
    local line
    while IFS= read -r line; do
        [[ -n $line ]] || continue
        out+=("$(shell_quote "$line")")
    done < "$file"
    printf '%s ' "${out[@]}"
}

time_field() {
    local key=$1 file=$2
    awk -F '\t' -v k="$key" '$1 == k { print $2; found = 1; exit } END { if (!found) print "NA" }' "$file"
}

run_case() {
    local case_id=$1 group=$2 engine=$3 scenario=$4 command=$5
    printf '%s\t%s\t%s\t%s\t%s\n' "$case_id" "$group" "$engine" "$scenario" "$command" >> "$COMMANDS_TSV"
    local repeat
    for repeat in $(seq 1 "$REPEATS"); do
        local stdout_file=$REPORT_DIR/${case_id}.r${repeat}.stdout
        local stderr_file=$REPORT_DIR/${case_id}.r${repeat}.stderr
        local time_file=$REPORT_DIR/${case_id}.r${repeat}.time
        local status=0
        log "run ${case_id} repeat ${repeat}/${REPEATS}"
        if [[ -n $TIME_BIN ]]; then
            "$TIME_BIN" -f 'elapsed_sec\t%e\nuser_sec\t%U\nsys_sec\t%S\nmax_rss_kb\t%M' \
                -o "$time_file" bash -lc "$command" > "$stdout_file" 2> "$stderr_file" || status=$?
        else
            local start_ns end_ns elapsed_ms
            start_ns=$(date +%s%N)
            bash -lc "$command" > "$stdout_file" 2> "$stderr_file" || status=$?
            end_ns=$(date +%s%N)
            elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
            {
                printf 'elapsed_sec\t%d.%03d\n' $((elapsed_ms / 1000)) $((elapsed_ms % 1000))
                printf 'user_sec\tNA\n'
                printf 'sys_sec\tNA\n'
                printf 'max_rss_kb\tNA\n'
            } > "$time_file"
        fi
        local elapsed user sys rss stdout_bytes stderr_bytes
        elapsed=$(time_field elapsed_sec "$time_file")
        user=$(time_field user_sec "$time_file")
        sys=$(time_field sys_sec "$time_file")
        rss=$(time_field max_rss_kb "$time_file")
        stdout_bytes=$(wc -c < "$stdout_file")
        stderr_bytes=$(wc -c < "$stderr_file")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$case_id" "$group" "$engine" "$scenario" "$repeat" \
            "$elapsed" "$user" "$sys" "$rss" "$stdout_bytes" "$stderr_bytes" \
            "$status" "$command" >> "$RAW_TSV"
        if [[ $status != 0 ]]; then
            echo "case failed: $case_id repeat $repeat" >&2
            sed -n '1,80p' "$stderr_file" >&2
            exit "$status"
        fi
    done
}

generate_synthetic_data() {
    log "generating synthetic benchmark data: N=$N"
    awk -v n="$N" 'BEGIN {
        print "@HD\tVN:1.6\tSO:coordinate";
        print "@SQ\tSN:chr1\tLN:" (n * 90 + 100000);
        print "@SQ\tSN:chr2\tLN:" (n * 90 + 100000);
        print "@SQ\tSN:chrX\tLN:" (n * 90 + 100000);
        print "@RG\tID:rg1\tLB:libA";
        print "@RG\tID:rg2\tLB:libB";
        seq = ""; qual = "";
        for (j = 0; j < 150; ++j) { seq = seq "A"; qual = qual "F"; }
        for (i = 1; i <= n; ++i) {
            chrom = (i % 3 == 0 ? "chrX" : (i % 2 == 0 ? "chr2" : "chr1"));
            rg = (i % 5 == 0 ? "rg2" : "rg1");
            pos = 1000 + i * 75;
            printf "q%08d\t0\t%s\t%d\t60\t150M\t*\t0\t0\t%s\t%s\tRG:Z:%s\tnp:i:%d\trq:f:0.99\n",
                   i, chrom, pos, seq, qual, rg, 1 + (i % 4);
        }
    }' > "$WORK/store.sam"
    awk -v n="$N" 'BEGIN {
        print "@HD\tVN:1.6\tSO:coordinate";
        print "@SQ\tSN:chr1\tLN:" (n * 90 + 100000);
        print "@SQ\tSN:chr2\tLN:" (n * 90 + 100000);
        print "@SQ\tSN:chrX\tLN:" (n * 90 + 100000);
        print "@RG\tID:rg1\tLB:libA";
        print "@RG\tID:rg2\tLB:libB";
        seq = ""; qual = "";
        for (j = 0; j < 150; ++j) { seq = seq "N"; qual = qual "!"; }
        for (i = 1; i <= n; ++i) {
            chrom = (i % 3 == 0 ? "chrX" : (i % 2 == 0 ? "chr2" : "chr1"));
            rg = (i % 5 == 0 ? "rg2" : "rg1");
            pos = 1200 + i * 75 + (i % 7);
            flag = 0;
            printf "q%08d\t%d\t%s\t%d\t40\t150M\t*\t0\t0\t%s\t%s\tRG:Z:%s\tNM:i:%d\tZZ:Z:%s\n",
                   i, flag, chrom, pos, seq, qual, rg, i % 3, (i % 11 == 0 ? "rare" : "common");
            if (i % 997 == 0) {
                printf "q%08d\t2048\t%s\t%d\t20\t75M75S\t*\t0\t0\t%s\t%s\tRG:Z:%s\tSA:Z:%s,%d,+,150M,40,0;\n",
                       i, chrom, pos + 20000, seq, qual, rg, chrom, pos;
            }
        }
    }' > "$WORK/aln.sam"
    "$SAMTOOLS" view -bS "$WORK/store.sam" -o "$WORK/store.unsorted.bam"
    "$SAMTOOLS" view -bS "$WORK/aln.sam" -o "$WORK/aln.unsorted.bam"
    "$SAMTOOLS" sort -o "$WORK/store.bam" "$WORK/store.unsorted.bam"
    "$SAMTOOLS" sort -o "$WORK/aln.bam" "$WORK/aln.unsorted.bam"
    "$SAMTOOLS" index "$WORK/store.bam"
    "$SAMTOOLS" index "$WORK/aln.bam"
    CANON_BAM=$WORK/store.bam
    ALN_BAM=$WORK/aln.bam
}

sample_real_data_if_requested() {
    if [[ $SAMPLE_N == 0 ]]; then
        return
    fi
    log "sampling $SAMPLE_N qnames from real canonical BAM"
    "$SAMTOOLS" view -F 2304 "$CANON_BAM" \
        | awk -v limit="$SAMPLE_N" '!seen[$1]++ { print $1; if (++n >= limit) exit }' \
        > "$WORK/sample.qnames"
    [[ -s $WORK/sample.qnames ]] || { echo "no sampled qnames" >&2; exit 1; }
    "$SAMTOOLS" view -b -N "$WORK/sample.qnames" -o "$WORK/store.sample.bam" "$CANON_BAM"
    "$SAMTOOLS" view -b -N "$WORK/sample.qnames" -o "$WORK/aln.sample.bam" "$ALN_BAM"
    "$SAMTOOLS" index "$WORK/store.sample.bam"
    "$SAMTOOLS" index "$WORK/aln.sample.bam"
    CANON_BAM=$WORK/store.sample.bam
    ALN_BAM=$WORK/aln.sample.bam
}

ensure_bam_index() {
    local bam=$1
    if [[ -f ${bam}.bai || -f ${bam%.bam}.bai || -f ${bam}.csi ]]; then
        return
    fi
    echo "missing BAM/CRAM index for $bam; create one or run with SAMPLE_N>0/SYNTHETIC=1" >&2
    exit 1
}

prepare_isam() {
    if [[ $PREPARE == 0 ]]; then
        [[ -n $ALN_ISAM && -f $ALN_ISAM ]] || { echo "PREPARE=0 requires ALN_ISAM" >&2; exit 1; }
        if [[ -z $STORE_SSAM ]]; then
            STORE_SSAM=${CANON_BAM%.*}.ssam
        fi
        return
    fi
    STORE_SSAM=$WORK/store.ssam
    ALN_ISAM=$WORK/aln.isam
    run_case prep_index prep isamtools index "$(shell_quote "$TOOL") index -o $(shell_quote "$STORE_SSAM") $(shell_quote "$CANON_BAM")"
    local import_cmd
    import_cmd="$(shell_quote "$TOOL") import $(shell_quote "$ALN_BAM") --store $(shell_quote "$CANON_BAM") --store-index $(shell_quote "$STORE_SSAM") --aln $(shell_quote "$ALN_ISAM")"
    if [[ $ALLOW_MISSING != 0 ]]; then
        import_cmd+=" --allow-missing"
    fi
    run_case prep_import prep isamtools import "$import_cmd"
    run_case prep_sort_all prep isamtools sort_all "$(shell_quote "$TOOL") sort --all $(shell_quote "$ALN_ISAM")"
}

make_query_inputs() {
    log "building qname, region, and BED inputs"
    "$SAMTOOLS" view -F 2304 "$ALN_BAM" \
        | awk -v limit="$QNAME_N" '!seen[$1]++ { print $1; if (++n >= limit) exit }' \
        > "$WORK/qnames.txt"
    [[ -s $WORK/qnames.txt ]] || { echo "no qnames selected from $ALN_BAM" >&2; exit 1; }

    "$SAMTOOLS" view -F 4 "$ALN_BAM" \
        | awk -v limit="$REGION_COUNT" -v width="$REGION_WINDOW" '
            $3 != "*" {
                beg = $4;
                end = beg + width;
                key = $3 ":" int((beg - 1) / width);
                if (!seen[key]++) {
                    print $3 ":" beg "-" end;
                    if (++n >= limit) exit;
                }
            }' > "$WORK/regions.txt"
    [[ -s $WORK/regions.txt ]] || { echo "no mapped regions selected from $ALN_BAM" >&2; exit 1; }
    head -n 1 "$WORK/regions.txt" > "$WORK/single_region.txt"
    awk -F '[:-]' 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $3 }' "$WORK/regions.txt" \
        | head -n "$TARGET_COUNT" > "$WORK/targets.bed"
    SINGLE_REGION=$(head -n 1 "$WORK/single_region.txt")
    MULTI_REGIONS=$(join_quoted_file "$WORK/regions.txt")
}

write_environment() {
    {
        echo "date_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "root $ROOT"
        echo "work $WORK"
        echo "tool $TOOL"
        "$TOOL" --version 2>/dev/null || true
        "$SAMTOOLS" --version 2>/dev/null | sed -n '1,4p' || true
        uname -a || true
        command -v lscpu >/dev/null && lscpu || true
        echo "synthetic $SYNTHETIC"
        echo "prepare $PREPARE"
        echo "repeats $REPEATS"
        echo "threads $THREADS"
        echo "time_bin ${TIME_BIN:-date-fallback}"
        echo "canonical_bam $CANON_BAM"
        echo "alignment_bam $ALN_BAM"
        echo "store_ssam ${STORE_SSAM:-NA}"
        echo "alignment_isam ${ALN_ISAM:-NA}"
    } > "$ENV_TXT"
}

collect_storage_stats() {
    {
        echo -e "path\tbytes"
        for path in "$CANON_BAM" "$ALN_BAM" "${STORE_SSAM:-}" "${ALN_ISAM:-}"; do
            [[ -n $path && -e $path ]] || continue
            printf '%s\t%s\n' "$path" "$(du -sb "$path" | awk '{ print $1 }')"
        done
    } > "$REPORT_DIR/file_sizes.tsv"
    if [[ -n ${STORE_SSAM:-} && -e $STORE_SSAM ]]; then
        "$TOOL" store_stats "$STORE_SSAM" > "$REPORT_DIR/store_stats.ssam.txt" || true
    fi
    if [[ -n ${ALN_ISAM:-} && -e $ALN_ISAM ]]; then
        "$TOOL" store_stats "$ALN_ISAM" > "$REPORT_DIR/store_stats.isam.txt" || true
    fi
    "$SAMTOOLS" idxstats "$ALN_BAM" > "$REPORT_DIR/alignment.idxstats.tsv" || true
}

compile_hts_benchmarks() {
    [[ $RUN_C_API != 0 ]] || return 0
    local htslib_prefix=${HTSLIB_PREFIX:-$ROOT/build/vendor/htslib}
    local lmdb_prefix=${LMDB_PREFIX:-$ROOT/build/vendor/lmdb}
    local bzip2_lib=${BZIP2_LIB:-$ROOT/build/vendor/bzip2/lib/libbz2.a}
    local htslib_inc=$htslib_prefix/include
    local htslib_lib=$htslib_prefix/lib/libhts.a
    local lmdb_lib=$lmdb_prefix/lib/liblmdb.a
    if [[ ! -f $htslib_lib || ! -f $lmdb_lib ]]; then
        echo "missing static htslib/lmdb build products; skipping C API benchmarks" >&2
        RUN_C_API=0
        return 0
    fi
    local bzip_link=()
    if [[ -f $bzip2_lib ]]; then
        bzip_link=("$bzip2_lib")
    else
        bzip_link=(-lbz2)
    fi
    cat > "$WORK/bench_hts.c" <<'C'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

static void die_msg(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(2);
}

static void count_record(bam1_t *b, long long *records, long long *bases) {
    ++*records;
    *bases += b->core.l_qseq;
}

static void run_full(const char *path) {
    htsFile *fp = sam_open(path, "r");
    if (!fp) die_msg("sam_open failed");
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) die_msg("sam_hdr_read failed");
    bam1_t *b = bam_init1();
    long long records = 0, bases = 0;
    int rc;
    while ((rc = sam_read1(fp, hdr, b)) >= 0)
        count_record(b, &records, &bases);
    if (rc < -1) die_msg("sam_read1 error");
    printf("records\t%lld\tbases\t%lld\n", records, bases);
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

static void run_regions(const char *path, const char *regions_path) {
    htsFile *fp = sam_open(path, "r");
    if (!fp) die_msg("sam_open failed");
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) die_msg("sam_hdr_read failed");
    hts_idx_t *idx = sam_index_load(fp, path);
    if (!idx) die_msg("sam_index_load failed");
    FILE *regions = fopen(regions_path, "r");
    if (!regions) die_msg("open regions failed");
    bam1_t *b = bam_init1();
    char *line = NULL;
    size_t cap = 0;
    long long records = 0, bases = 0, region_count = 0;
    while (getline(&line, &cap, regions) >= 0) {
        line[strcspn(line, "\r\n")] = 0;
        if (!line[0]) continue;
        hts_itr_t *itr = sam_itr_querys(idx, hdr, line);
        if (!itr) die_msg("sam_itr_querys failed");
        int rc;
        while ((rc = sam_itr_next(fp, itr, b)) >= 0)
            count_record(b, &records, &bases);
        if (rc < -1) die_msg("sam_itr_next error");
        hts_itr_destroy(itr);
        ++region_count;
    }
    printf("regions\t%lld\trecords\t%lld\tbases\t%lld\n", region_count, records, bases);
    free(line);
    fclose(regions);
    bam_destroy1(b);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: bench_hts full <bam-or-isam> | regions <bam-or-isam> <regions.txt>\n");
        return 2;
    }
    if (strcmp(argv[1], "full") == 0) {
        run_full(argv[2]);
    } else if (strcmp(argv[1], "regions") == 0) {
        if (argc != 4) return 2;
        run_regions(argv[2], argv[3]);
    } else {
        return 2;
    }
    return 0;
}
C
    "$CC" -O3 -g -no-pie -I"$htslib_inc" "$WORK/bench_hts.c" -o "$WORK/bench_hts_native" \
        "$htslib_lib" "${bzip_link[@]}" -lz -llzma -lm -lpthread -lzstd \
        2> "$REPORT_DIR/compile.htslib-native.stderr" || {
            cat "$REPORT_DIR/compile.htslib-native.stderr" >&2
            exit 1
        }
    "$CC" -O3 -g -no-pie -I"$ROOT/include" -I"$htslib_inc" -include htsisam/sam.h \
        "$WORK/bench_hts.c" -o "$WORK/bench_hts_compat" \
        "$ROOT/build/libisam.a" "$lmdb_lib" "$htslib_lib" "${bzip_link[@]}" \
        -lz -llzma -lm -lpthread -lzstd \
        2> "$REPORT_DIR/compile.htsisam.stderr" || {
            cat "$REPORT_DIR/compile.htsisam.stderr" >&2
            exit 1
        }
}

run_cli_benchmarks() {
    local st="$SAMTOOLS"
    local it="$TOOL"
    local bam
    local isam
    bam=$(shell_quote "$ALN_BAM")
    isam=$(shell_quote "$ALN_ISAM")
    local qnames targets single
    qnames=$(shell_quote "$WORK/qnames.txt")
    targets=$(shell_quote "$WORK/targets.bed")
    single=$(shell_quote "$SINGLE_REGION")

    run_case cli_samtools_count_full cli_view samtools count_full \
        "$(shell_quote "$st") view -@ $THREADS -c $bam"
    run_case cli_isamtools_count_full cli_view isamtools count_full \
        "$(shell_quote "$it") view -c $isam"

    run_case cli_samtools_sam_full cli_view samtools sam_full \
        "$(shell_quote "$st") view -@ $THREADS $bam > /dev/null"
    run_case cli_isamtools_sam_full cli_view isamtools sam_full \
        "$(shell_quote "$it") view --no-header $isam > /dev/null"

    run_case cli_samtools_count_region cli_view samtools count_single_region \
        "$(shell_quote "$st") view -@ $THREADS -c $bam $single"
    run_case cli_isamtools_count_region cli_view isamtools count_single_region \
        "$(shell_quote "$it") view -c $isam $single"

    run_case cli_samtools_count_multiregion cli_view samtools count_multi_region \
        "$(shell_quote "$st") view -@ $THREADS -c $bam $MULTI_REGIONS"
    run_case cli_isamtools_count_multiregion cli_view isamtools count_multi_region \
        "$(shell_quote "$it") view -c $isam $MULTI_REGIONS"

    run_case cli_samtools_count_bed cli_view samtools count_bed_targets \
        "$(shell_quote "$st") view -@ $THREADS -c -L $targets $bam"
    run_case cli_isamtools_count_bed cli_view isamtools count_bed_targets \
        "$(shell_quote "$it") view -c -L $targets $isam"

    run_case cli_samtools_count_qnames cli_view samtools count_qnames \
        "$(shell_quote "$st") view -@ $THREADS -c -N $qnames $bam"
    run_case cli_isamtools_count_qnames cli_view isamtools count_qnames \
        "$(shell_quote "$it") view -c -N $qnames $isam"

    run_case cli_samtools_count_filter cli_view samtools count_filter_F2304_q20 \
        "$(shell_quote "$st") view -@ $THREADS -c -F 2304 -q 20 $bam"
    run_case cli_isamtools_count_filter cli_view isamtools count_filter_F2304_q20 \
        "$(shell_quote "$it") view -c -F 2304 -q 20 $isam"

    run_case cli_samtools_bam_region cli_view samtools bam_single_region \
        "$(shell_quote "$st") view -@ $THREADS -b -o /dev/null $bam $single"
    run_case cli_isamtools_bam_region cli_view isamtools bam_single_region \
        "$(shell_quote "$it") view -b -o /dev/null $isam $single"

    run_case cli_samtools_fast_bam_region cli_view samtools fast_bam_single_region \
        "$(shell_quote "$st") view -@ $THREADS -1 -o /dev/null $bam $single"
    run_case cli_isamtools_fast_bam_region cli_view isamtools fast_bam_single_region \
        "$(shell_quote "$it") view -1 -o /dev/null $isam $single"

    if [[ -n ${STORE_SSAM:-} && -e $STORE_SSAM ]]; then
        local ssam region
        ssam=$(shell_quote "$STORE_SSAM")
        region=$(shell_quote "$SINGLE_REGION")
        run_case query_region query isamtools query_region \
            "$(shell_quote "$it") query --source store:$ssam --source aln:$isam --select aln:$region --columns qname,store.region,aln.region > /dev/null"
        run_case query_span query isamtools query_span \
            "$(shell_quote "$it") query --source store:$ssam --source aln:$isam --select aln:$region --span aln:$SPAN_DIST --columns qname,store.region,aln.region,aln.n > /dev/null"
    fi
}

run_c_api_benchmarks() {
    [[ $RUN_C_API != 0 ]] || return 0
    local bam isam regions
    bam=$(shell_quote "$ALN_BAM")
    isam=$(shell_quote "$ALN_ISAM")
    regions=$(shell_quote "$WORK/regions.txt")
    run_case c_native_bam_full c_api htslib-native full_bam \
        "$(shell_quote "$WORK/bench_hts_native") full $bam"
    run_case c_compat_bam_full c_api htsisam-dispatch full_bam \
        "$(shell_quote "$WORK/bench_hts_compat") full $bam"
    run_case c_compat_isam_full c_api htsisam-isam full_isam \
        "$(shell_quote "$WORK/bench_hts_compat") full $isam"
    run_case c_native_bam_regions c_api htslib-native multi_region_bam \
        "$(shell_quote "$WORK/bench_hts_native") regions $bam $regions"
    run_case c_compat_bam_regions c_api htsisam-dispatch multi_region_bam \
        "$(shell_quote "$WORK/bench_hts_compat") regions $bam $regions"
    run_case c_compat_isam_regions c_api htsisam-isam multi_region_isam \
        "$(shell_quote "$WORK/bench_hts_compat") regions $isam $regions"
}

run_perf_samples() {
    [[ $RUN_PERF != 0 ]] || return 0
    command -v perf >/dev/null || { echo "perf not found; skipping RUN_PERF=1 samples" >&2; return 0; }
    local bam isam single
    bam=$(shell_quote "$ALN_BAM")
    isam=$(shell_quote "$ALN_ISAM")
    single=$(shell_quote "$SINGLE_REGION")
    perf stat -d -o "$REPORT_DIR/perf.samtools_region.txt" \
        "$SAMTOOLS" view -c "$ALN_BAM" "$SINGLE_REGION" >/dev/null || true
    perf stat -d -o "$REPORT_DIR/perf.isamtools_region.txt" \
        "$TOOL" view -c "$ALN_ISAM" "$SINGLE_REGION" >/dev/null || true
    if [[ $RUN_C_API != 0 ]]; then
        bash -lc "perf stat -d -o $(shell_quote "$REPORT_DIR/perf.htsisam_isam_regions.txt") $(shell_quote "$WORK/bench_hts_compat") regions $isam $(shell_quote "$WORK/regions.txt") >/dev/null" || true
        bash -lc "perf stat -d -o $(shell_quote "$REPORT_DIR/perf.htslib_native_regions.txt") $(shell_quote "$WORK/bench_hts_native") regions $bam $(shell_quote "$WORK/regions.txt") >/dev/null" || true
    fi
}

write_summary() {
    local tmp=$REPORT_DIR/summary.unsorted.tsv
    awk -F '\t' '
        BEGIN {
            OFS = "\t";
            print "case_id", "group", "engine", "scenario", "n", "mean_ms", "min_ms", "mean_user_ms", "mean_sys_ms", "mean_max_rss_kb";
        }
        NR == 1 { next }
        $12 != 0 || $6 == "NA" { next }
        {
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4;
            n[key]++;
            elapsed[key] += $6 * 1000.0;
            user[key] += $7 * 1000.0;
            sys[key] += $8 * 1000.0;
            rss[key] += $9;
            if (!(key in min) || $6 * 1000.0 < min[key]) min[key] = $6 * 1000.0;
            case_id[key] = $1;
            group[key] = $2;
            engine[key] = $3;
            scenario[key] = $4;
        }
        END {
            for (key in n) {
                print case_id[key], group[key], engine[key], scenario[key], n[key],
                      sprintf("%.3f", elapsed[key] / n[key]),
                      sprintf("%.3f", min[key]),
                      sprintf("%.3f", user[key] / n[key]),
                      sprintf("%.3f", sys[key] / n[key]),
                      sprintf("%.0f", rss[key] / n[key]);
            }
        }' "$RAW_TSV" > "$tmp"
    {
        head -n 1 "$tmp"
        tail -n +2 "$tmp" | sort -t $'\t' -k2,2 -k4,4 -k3,3
    } > "$SUMMARY_TSV"
}

write_report() {
    {
        echo "# isamtools Performance Profile"
        echo
        echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo
        echo "## Benchmark Design"
        echo
        echo "- CLI comparisons use \`samtools view\` on BAM versus \`isamtools view\` on the imported alignment view."
        echo "- C API comparisons use a tiny htslib-style reader compiled two ways: native htslib, and the htsisam drop-in dispatch layer."
        echo "- \`htsisam-dispatch\` on BAM estimates the wrapper overhead on ordinary BAM input."
        echo "- \`htsisam-isam\` on .isam estimates the replacement path seen by htslib consumers."
        echo "- Timings are wall-clock means over \`REPEATS=$REPEATS\`; raw per-repeat timings are in \`raw.tsv\`."
        echo "- If this was run over a network mount, treat absolute wall times as noisy; ratios across adjacent cases are usually more useful."
        echo
        echo "## Inputs"
        echo
        echo '```'
        sed -n '1,80p' "$ENV_TXT"
        echo '```'
        echo
        echo "## File Sizes"
        echo
        echo '| path | bytes |'
        echo '|---|---:|'
        awk -F '\t' 'NR > 1 { printf "| `%s` | %s |\n", $1, $2 }' "$REPORT_DIR/file_sizes.tsv"
        echo
        echo "## Summary"
        echo
        echo '| group | scenario | engine | n | mean ms | min ms | user ms | sys ms | RSS KB | vs baseline |'
        echo '|---|---|---|---:|---:|---:|---:|---:|---:|---:|'
        {
        awk -F '\t' '
            function baseline_key(group, scenario, s) {
                s = scenario;
                if (group == "c_api")
                    sub(/_isam$/, "_bam", s);
                return group SUBSEP s;
            }
            NR == 1 { next }
            {
                rows[++nr] = $0;
                key = $2 SUBSEP $4;
                if (($3 == "samtools" || $3 == "htslib-native") && !(key in base))
                    base[key] = $6;
            }
            END {
                for (i = 1; i <= nr; ++i) {
                    split(rows[i], f, FS);
                    key = baseline_key(f[2], f[4]);
                    ratio = (key in base && base[key] > 0) ? f[6] / base[key] : 1.0;
                    printf "| %s | %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %.0f | %.2fx |\n",
                           f[2], f[4], f[3], f[5], f[6], f[7], f[8], f[9], f[10], ratio;
                }
            }' "$SUMMARY_TSV"
        echo
        echo "## Largest Relative Slowdowns"
        echo
        echo '| group | scenario | engine | mean ms | baseline ms | ratio |'
        echo '|---|---|---|---:|---:|---:|'
        awk -F '\t' '
            function baseline_key(group, scenario, s) {
                s = scenario;
                if (group == "c_api")
                    sub(/_isam$/, "_bam", s);
                return group SUBSEP s;
            }
            NR == 1 { next }
            {
                rows[++nr] = $0;
                key = $2 SUBSEP $4;
                if (($3 == "samtools" || $3 == "htslib-native") && !(key in base))
                    base[key] = $6;
            }
            END {
                for (i = 1; i <= nr; ++i) {
                    split(rows[i], f, FS);
                    key = baseline_key(f[2], f[4]);
                    if (!(key in base) || base[key] <= 0 || f[3] == "samtools" || f[3] == "htslib-native")
                        continue;
                    ratio = f[6] / base[key];
                    printf "%.6f\t| %s | %s | %s | %.3f | %.3f | %.2fx |\n",
                           ratio, f[2], f[4], f[3], f[6], base[key], ratio;
                }
            }' "$SUMMARY_TSV" | sort -nr | head -n 12 | cut -f2- || true
        }
        echo
        echo "## Bottleneck Reading Guide"
        echo
        echo "- \`count_full\` mostly measures iteration, filtering, and record reconstruction without formatting output."
        echo "- \`sam_full\` adds SAM formatting cost. A large gap from \`count_full\` points at formatting or reconstruction of emitted fields."
        echo "- \`bam_single_region\` and \`fast_bam_single_region\` add BAM encoding/compression. If \`fast_bam\` is much faster, compression is material."
        echo "- \`count_single_region\`, \`count_multi_region\`, and \`count_bed_targets\` stress interval lookup and repeated iterator creation."
        echo "- \`count_qnames\` stresses name filtering; for isamtools this should expose qname-index and store reconstruction costs."
        echo "- \`htsisam-dispatch\` versus \`htslib-native\` on BAM is the pure drop-in wrapper overhead."
        echo "- \`htsisam-isam\` versus \`htslib-native\` is the end-to-end cost of replacing BAM input with .isam for an htslib consumer."
        echo "- High mean system time relative to user time usually points at IO, LMDB paging, or network filesystem effects."
        echo
        echo "## Artifacts"
        echo
        echo "- Raw timings: \`$RAW_TSV\`"
        echo "- Summary table: \`$SUMMARY_TSV\`"
        echo "- Commands: \`$COMMANDS_TSV\`"
        echo "- Environment: \`$ENV_TXT\`"
        if [[ -f $REPORT_DIR/store_stats.isam.txt ]]; then
            echo "- Alignment store stats: \`$REPORT_DIR/store_stats.isam.txt\`"
        fi
        if [[ -f $REPORT_DIR/store_stats.ssam.txt ]]; then
            echo "- Store stats: \`$REPORT_DIR/store_stats.ssam.txt\`"
        fi
        if [[ $RUN_PERF != 0 ]]; then
            echo "- Perf samples: \`$REPORT_DIR/perf.*.txt\`"
        fi
    } > "$REPORT_MD"
}

if [[ $SYNTHETIC != 0 ]]; then
    generate_synthetic_data
else
    [[ -f $CANON_BAM ]] || { echo "missing canonical BAM: $CANON_BAM" >&2; exit 1; }
    [[ -f $ALN_BAM ]] || { echo "missing alignment BAM: $ALN_BAM" >&2; exit 1; }
    sample_real_data_if_requested
    ensure_bam_index "$CANON_BAM"
    ensure_bam_index "$ALN_BAM"
fi

prepare_isam
write_environment
make_query_inputs
collect_storage_stats
compile_hts_benchmarks
run_cli_benchmarks
run_c_api_benchmarks
run_perf_samples
write_summary
write_report

log "raw timings: $RAW_TSV"
log "summary: $SUMMARY_TSV"
log "report: $REPORT_MD"
