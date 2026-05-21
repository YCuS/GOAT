#!/usr/bin/env bash
set -euo pipefail

cfg="${1:-config.json}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

command -v jq >/dev/null 2>&1 || {
  echo "ERROR: jq is required to read config.json" >&2
  exit 1
}

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

orig_pwm=$(jq -r '.paths.original_pwm_file' "$cfg")
seq_path=$(jq -r '.paths.sequence_file' "$cfg")
res_path=$(jq -r '.paths.results_file' "$cfg")

bayes_raw=$(jq -r '.parameters.bayes // false' "$cfg")
if [[ "$bayes_raw" == "true" || "$bayes_raw" == "1" ]]; then
  bayes=1
else
  bayes=0
fi

pseudocount=$(jq -r '.parameters.pwm2fm // 45' "$cfg")
full_pct=$(jq -r '.parameters.thr1 // 95.0' "$cfg")
core_pct=$(jq -r '.parameters.thr2 // 75.0' "$cfg")
iterations=$(jq -r '.parameters.simulation_iterations // 10000' "$cfg")
core_len=$(jq -r '.settings.core_length // 5' "$cfg")
precision=$(jq -r '.settings.precision // 6' "$cfg")

save_thr_raw=$(jq -r '.settings.save_threshold // true' "$cfg")
threshold_from_cfg=$(jq -r '.paths.threshold_file // empty' "$cfg")

if [[ -n "$threshold_from_cfg" ]]; then
  threshold_out="$threshold_from_cfg"
elif [[ "$save_thr_raw" == "true" || "$save_thr_raw" == "1" ]]; then
  if [[ "$res_path" == *.txt || "$res_path" == *.tsv ]]; then
    threshold_out="$(dirname -- "$res_path")/thresholds.tsv"
  else
    threshold_out="${res_path%/}/thresholds.tsv"
  fi
else
  threshold_out="$tmp_dir/thresholds.tsv"
fi

pwm_out="$tmp_dir/pwm.tsv"
mkdir -p "$(dirname -- "$threshold_out")"

echo "[1/3] Converting motifs to normalized PWM"
./meme2pwm "$orig_pwm" "$pwm_out" "$bayes" "$pseudocount"

echo "[2/3] Calculating fixed PWM thresholds"
./get_thr "$pwm_out" "$threshold_out" "$full_pct" "$core_pct" "$iterations" "$core_len" "$precision"

echo "[3/3] Searching motif hits"
if [[ -d "$seq_path" ]]; then
  if [[ "$res_path" == *.txt || "$res_path" == *.tsv ]]; then
    output_path="${res_path%.*}"
  else
    output_path="${res_path%/}"
  fi
  mkdir -p "$output_path"
  ./search_motif "$seq_path" "$pwm_out" "$threshold_out" "$output_path" "$core_len" --batch
  echo "Results directory: $output_path"
else
  mkdir -p "$(dirname -- "$res_path")"
  ./search_motif "$seq_path" "$pwm_out" "$threshold_out" "$res_path" "$core_len"
  echo "Results file: $res_path"
fi

if [[ "$threshold_out" == "$tmp_dir/"* ]]; then
  echo "Thresholds: temporary only"
else
  echo "Thresholds: $threshold_out"
fi
echo "Done."
