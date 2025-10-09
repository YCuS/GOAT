#!/usr/bin/env bash
set -euo pipefail

cfg="config.json"

# read paths from config.json
orig_pwm=$(jq -r .paths.original_pwm_file "$cfg")
seq_path=$(jq -r .paths.sequence_file     "$cfg")  # may be FASTA file or directory
res_path=$(jq -r .paths.results_file       "$cfg")  # file or directory
allow_gaps_raw=$(jq -r '.settings.allow_gaps // "true"' "$cfg")
if [[ "$allow_gaps_raw" == "true" || "$allow_gaps_raw" == "1" ]]; then
  gaps_flag=""
else
  gaps_flag="--no-gaps"
fi

# read settings: whether to persist thresholds (default: true)
save_thr_raw=$(jq -r '.settings.save_threshold // "true"' "$cfg")
if [[ "$save_thr_raw" == "true" || "$save_thr_raw" == "1" ]]; then
  save_thr=1
else
  save_thr=0
fi

# create temporary working directory early (may be used for temp thr)
TMP=$(mktemp -d)

# resolve threshold output path
# precedence: explicit .paths.threshold_file (force save) > save_thr flag > TMP
thr_path_from_cfg=$(jq -r '.paths.threshold_file // empty' "$cfg")
if [[ -n "${thr_path_from_cfg:-}" ]]; then
  thr_out="$thr_path_from_cfg"        # explicit path -> always persist
  mkdir -p "$(dirname -- "$thr_out")"
elif [[ $save_thr -eq 1 ]]; then
  if [[ "$res_path" == *.txt ]]; then
    thr_dir=$(dirname -- "$res_path")
    thr_out="${thr_dir}/thr.txt"
  else
    thr_out="${res_path%/}/thr.txt"
  fi
  mkdir -p "$(dirname -- "$thr_out")"
else
  thr_out="$TMP/thr.txt"              # no persist -> temp file in TMP
fi

# read parameters
bayes_raw=$(jq -r .parameters.bayes "$cfg")
if [[ "$bayes_raw" == "true" ]]; then bayes=1; else bayes=0; fi
count=$(jq -r .parameters.pwm2fm               "$cfg")
thr1=$(jq -r .parameters.thr1                  "$cfg")
thr2=$(jq -r .parameters.thr2                  "$cfg")
iters=$(jq -r .parameters.simulation_iterations "$cfg")
core=$(jq -r .settings.core_length             "$cfg")
prec=$(jq -r .settings.precision               "$cfg")

# Optional advanced threshold controls
gc_thr=$(jq -r '.settings.gc_threshold // empty' "$cfg")
rel_full=$(jq -r '.parameters.relaxed_full // empty' "$cfg")
rel_core=$(jq -r '.parameters.relaxed_core // empty' "$cfg")
str_full=$(jq -r '.parameters.strict_full // empty' "$cfg")
str_core=$(jq -r '.parameters.strict_core // empty' "$cfg")
cluster_enable=$(jq -r '.clustering.enable // empty' "$cfg")
cluster_k=$(jq -r '.clustering.k // empty' "$cfg")
cluster_relaxed=$(jq -r '.clustering.relaxed_cluster // empty' "$cfg")

# 1) convert MEME motifs to PWM
./meme2pwm "$orig_pwm" "$TMP/pwm.txt" "$bayes" "$count"

# 2) calculate thresholds -> write to thr_out (persist or temp)
GET_THR_FLAGS=()
[[ -n "$gc_thr"   && "$gc_thr"   != "null" ]] && GET_THR_FLAGS+=("--gc-thresh=$gc_thr")
[[ -n "$rel_full" && "$rel_full" != "null" ]] && GET_THR_FLAGS+=("--relaxed-full=$rel_full")
[[ -n "$rel_core" && "$rel_core" != "null" ]] && GET_THR_FLAGS+=("--relaxed-core=$rel_core")
[[ -n "$str_full" && "$str_full" != "null" ]] && GET_THR_FLAGS+=("--strict-full=$str_full")
[[ -n "$str_core" && "$str_core" != "null" ]] && GET_THR_FLAGS+=("--strict-core=$str_core")

./get_thr "$TMP/pwm.txt" "$thr_out" "$thr1" "$thr2" "$iters" "$core" "$prec" "${GET_THR_FLAGS[@]}"
# Append clustering flags if set
EXTRA_CLUSTER_FLAGS=()
if [[ -n "$cluster_enable" && "$cluster_enable" != "null" ]]; then
  EXTRA_CLUSTER_FLAGS+=("--cluster-enable=$cluster_enable")
fi
if [[ -n "$cluster_k" && "$cluster_k" != "null" ]]; then
  EXTRA_CLUSTER_FLAGS+=("--cluster-k=$cluster_k")
fi
if [[ -n "$cluster_relaxed" && "$cluster_relaxed" != "null" ]]; then
  EXTRA_CLUSTER_FLAGS+=("--cluster-relaxed=$cluster_relaxed")
fi

if (( ${#EXTRA_CLUSTER_FLAGS[@]} )); then
  # derive plot path (same directory as thr_out)
  plot_dir=$(dirname -- "$thr_out")
  plot_path="$plot_dir/cluster.svg"
  EXTRA_CLUSTER_FLAGS+=("--cluster-plot=$plot_path")
  ./get_thr "$TMP/pwm.txt" "$thr_out" "$thr1" "$thr2" "$iters" "$core" "$prec" "${GET_THR_FLAGS[@]}" "${EXTRA_CLUSTER_FLAGS[@]}" >/dev/null 2>&1 || true
  echo "Cluster plot (if clustering enabled) saved to: $plot_path"
fi

# verify intermediate files
[[ -s "$TMP/pwm.txt" ]] || { echo "ERROR: pwm.txt is empty"; exit 1; }
[[ -s "$thr_out"    ]] || { echo "ERROR: thr.txt is empty ($thr_out)"; exit 1; }

# 3) search motifs
if [[ -d "$seq_path" ]]; then
  if [[ "$res_path" == *.txt ]]; then out_dir="${res_path%.*}"; else out_dir="$res_path"; fi
  mkdir -p "$out_dir"
  ./search_motif "$seq_path" "$TMP/pwm.txt" "$thr_out" "$out_dir" "$core" --batch ${gaps_flag:+$gaps_flag}
  if [[ -n "${thr_path_from_cfg:-}" || $save_thr -eq 1 ]]; then
    echo "Batch search completed. Thresholds saved to: $thr_out"
  else
    echo "Batch search completed. Thresholds were temporary (not saved)."
  fi
  echo "Results saved under directory: $out_dir"
else
  ./search_motif "$seq_path" "$TMP/pwm.txt" "$thr_out" "$res_path" "$core" ${gaps_flag:+$gaps_flag}
  if [[ -n "${thr_path_from_cfg:-}" || $save_thr -eq 1 ]]; then
    echo "Search completed. Thresholds saved to: $thr_out"
  else
    echo "Search completed. Thresholds were temporary (not saved)."
  fi
  echo "Results saved to file: $res_path"
fi

# clean up (TMP removal will delete temp thr.txt if not persisted)
rm -rf "$TMP"