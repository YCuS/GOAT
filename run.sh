#!/usr/bin/env bash
set -euo pipefail

cfg="config.json"

# read paths
orig_pwm=$(jq -r .paths.original_pwm_file "$cfg")
seq_fasta=$(jq -r .paths.sequence_file     "$cfg")
res_file=$(jq -r .paths.results_file       "$cfg")

# read parameters
bayes_raw=$(jq -r .parameters.bayes "$cfg")
# convert JSON true/false → 1/0
if [[ "$bayes_raw" == "true" ]]; then
  bayes=1
else
  bayes=0
fi

count=$(jq -r .parameters.pwm2fm "$cfg")
thr1=$(jq -r .parameters.thr1  "$cfg")
thr2=$(jq -r .parameters.thr2  "$cfg")
iters=$(jq -r .parameters.simulation_iterations "$cfg")
core=$(jq -r .settings.core_length "$cfg")
prec=$(jq -r .settings.precision   "$cfg")

# make temp dir
TMP=$(mktemp -d)

# 1) MEME → PWM
./meme2pwm "$orig_pwm" "$TMP/pwm.txt" "$bayes" "$count"

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