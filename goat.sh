#!/usr/bin/env bash
set -euo pipefail

# 默认值
thr1=95
thr2=75
num=45
core=5
bayes_raw=true
iters=10000
prec=6

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -pwm)
            orig_pwm="$2"
            shift 2
            ;;
        -seq)
            seq_fasta="$2"
            shift 2
            ;;
        -out)
            res_file="$2"
            shift 2
            ;;
        -thr1)
            thr1="$2"
            shift 2
            ;;
        -thr2)
            thr2="$2"
            shift 2
            ;;
        -pwm2fm)
            num="$2"
            shift 2
            ;;
        -core)
            core="$2"
            shift 2
            ;;
        -bayes)
            bayes_raw="$2"
            shift 2
            ;;
        -sim)
            iters="$2"
            shift 2
            ;;
        -precsie)
            prec="$2"
            shift 2
            ;;
        *)
            echo "未知选项: $1" >&2
            exit 1
            ;;
    esac
done

# 检查必需参数
if [ -z "${orig_pwm:-}" ] || [ -z "${seq_fasta:-}" ] || [ -z "${res_file:-}" ]; then
    echo "error: -pwm、-seq and -out option is necessary" >&2
    exit 1
fi

# convert JSON true/false → 1/0
if [[ "$bayes_raw" == "true" ]]; then
    bayes=1
else
    bayes=0
fi

# make temp dir
TMP=$(mktemp -d)

# 1) MEME → PWM
./meme2pwm "$orig_pwm" "$TMP/pwm.txt" "$bayes" "$num"

# 2) calculate thresholds
./get_thr "$TMP/pwm.txt" "$TMP/thr.txt" \
          "$thr1" "$thr2" "$iters" "$core" "$prec"

# check intermediate files
[[ -s "$TMP/pwm.txt" ]] || { echo "ERROR: pwm.txt is empty"; exit 1; }
[[ -s "$TMP/thr.txt" ]] || { echo "ERROR: thr.txt is empty"; exit 1; }

# 3) search motifs
./search_motif "$seq_fasta" "$TMP/pwm.txt" "$TMP/thr.txt" "$res_file" "$core"

echo "Done. Results → $res_file"

# clean up
rm -rf "$TMP"
