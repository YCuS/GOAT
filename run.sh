#!/usr/bin/env bash
set -euo pipefail

cfg="config.json"

# read paths from config.json
orig_pwm=$(jq -r .paths.original_pwm_file "$cfg")
seq_path=$(jq -r .paths.sequence_file     "$cfg")  # may be FASTA file or directory
res_path=$(jq -r .paths.results_file       "$cfg")  # file or directory

# read parameters
bayes_raw=$(jq -r .parameters.bayes "$cfg")
if [[ "$bayes_raw" == "true" ]]; then
  bayes=1
else
  bayes=0
fi

count=$(jq -r .parameters.pwm2fm               "$cfg")
thr1=$(jq -r .parameters.thr1                  "$cfg")
thr2=$(jq -r .parameters.thr2                  "$cfg")
iters=$(jq -r .parameters.simulation_iterations "$cfg")
core=$(jq -r .settings.core_length             "$cfg")
prec=$(jq -r .settings.precision               "$cfg")

# create temporary working directory
TMP=$(mktemp -d)

# 1) convert MEME motifs to PWM
./meme2pwm "$orig_pwm" "$TMP/pwm.txt" "$bayes" "$count"

# 2) calculate thresholds
./get_thr "$TMP/pwm.txt" "$TMP/thr.txt" \
          "$thr1" "$thr2" "$iters" "$core" "$prec"

# verify intermediate files
[[ -s "$TMP/pwm.txt" ]] || { echo "ERROR: pwm.txt is empty"; exit 1; }
[[ -s "$TMP/thr.txt" ]] || { echo "ERROR: thr.txt is empty"; exit 1; }

# 3) search motifs: batch if seq_path is directory, otherwise single-file mode
if [[ -d "$seq_path" ]]; then
  # batch mode: decide output directory
  if [[ "$res_path" == *.txt ]]; then
    # strip extension to create a directory
    out_dir="${res_path%.*}"
  else
    out_dir="$res_path"
  fi
  mkdir -p "$out_dir"
  ./search_motif "$seq_path" "$TMP/pwm.txt" "$TMP/thr.txt" "$out_dir" "$core" --batch
  echo "Batch search completed. Results saved under directory: $out_dir"
else
  # single-file mode: results path is a file
  ./search_motif "$seq_path" "$TMP/pwm.txt" "$TMP/thr.txt" "$res_path" "$core"
  echo "Search completed. Results saved to file: $res_path"
fi

# clean up
rm -rf "$TMP"