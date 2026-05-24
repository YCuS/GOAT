#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash run.sh <motif_file> <sequence_file_or_dir> <results_file_or_dir> [options]

Required arguments:
  motif_file           MEME file or simple Motif: file
  sequence_file_or_dir FASTA file, or directory with .fa/.fasta/.fna files
  results_file_or_dir  output TSV file; output directory when sequence input is a directory

Options:
  --threshold-file PATH  save thresholds here
  --core-length N        core window length (default: 5)
  --full-pct PCT         full motif threshold percentile (default: 95.0)
  --core-pct PCT         core threshold percentile (default: 75.0)
  --iterations N         simulations per motif (default: 10000)
  --precision N          threshold decimal precision (default: 6)
  --bayes 0|1            Bayesianize/smooth PWM probabilities in meme2pwm (default: 0)
  --pseudocount N        pseudocount value passed to meme2pwm (default: 45)
  -h, --help             show this help

Single-file example:
  bash run.sh examples/motifs.txt examples/sequences.fa results/motif_hits.tsv

Batch example:
  bash run.sh examples/motifs.txt path/to/fasta_dir results/batch_hits
USAGE
}

if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 3 ]]; then
  usage >&2
  exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

motif_file="$1"
seq_path="$2"
res_path="$3"
shift 3

threshold_out=""
core_len=5
full_pct=95.0
core_pct=75.0
iterations=10000
precision=6
bayes=0
pseudocount=45

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threshold-file)
      threshold_out="${2:?ERROR: --threshold-file requires a path}"
      shift 2
      ;;
    --core-length)
      core_len="${2:?ERROR: --core-length requires a value}"
      shift 2
      ;;
    --full-pct)
      full_pct="${2:?ERROR: --full-pct requires a value}"
      shift 2
      ;;
    --core-pct)
      core_pct="${2:?ERROR: --core-pct requires a value}"
      shift 2
      ;;
    --iterations)
      iterations="${2:?ERROR: --iterations requires a value}"
      shift 2
      ;;
    --precision)
      precision="${2:?ERROR: --precision requires a value}"
      shift 2
      ;;
    --bayes)
      bayes="${2:?ERROR: --bayes requires 0 or 1}"
      shift 2
      ;;
    --pseudocount)
      pseudocount="${2:?ERROR: --pseudocount requires a value}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[[ -f "$motif_file" ]] || { echo "ERROR: motif file not found: $motif_file" >&2; exit 1; }
[[ -e "$seq_path" ]] || { echo "ERROR: sequence path not found: $seq_path" >&2; exit 1; }

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

if [[ -z "$threshold_out" ]]; then
  if [[ -d "$seq_path" ]]; then
    threshold_out="${res_path%/}/thresholds.tsv"
  else
    threshold_out="$(dirname -- "$res_path")/thresholds.tsv"
  fi
fi

pwm_out="$tmp_dir/pwm.tsv"
mkdir -p "$(dirname -- "$threshold_out")"

echo "[1/3] Converting motifs to normalized PWM"
./meme2pwm "$motif_file" "$pwm_out" "$bayes" "$pseudocount"

echo "[2/3] Calculating fixed PWM thresholds"
./get_thr "$pwm_out" "$threshold_out" "$full_pct" "$core_pct" "$iterations" "$core_len" "$precision"

echo "[3/3] Searching motif hits"
if [[ -d "$seq_path" ]]; then
  mkdir -p "$res_path"
  ./search_motif "$seq_path" "$pwm_out" "$threshold_out" "$res_path" "$core_len" --batch
  echo "Results directory: $res_path"
else
  mkdir -p "$(dirname -- "$res_path")"
  ./search_motif "$seq_path" "$pwm_out" "$threshold_out" "$res_path" "$core_len"
  echo "Results file: $res_path"
fi

echo "Thresholds: $threshold_out"
echo "Done."
