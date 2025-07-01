# Goodness Of Alignment Test (GOAT)

This is a tool using GOAT (Goodness Of Alignment Test) to search for motifs in DNA(RNA) sequences using Position Weight Matrices (PWMs).

---

## Prerequisites

- macOS or Linux
- A C++17 compiler (e.g. `g++`)
- `jq` (for JSON parsing)
- Bash or Fish shell

---

## Repository Structure

```
compile.sh           # Build all binaries
run.sh               # Orchestrate end‐to‐end workflow
config.json          # User parameters and paths
meme2pwm.cpp         # Other format PWM→The PWM we use converter
get_thr.cpp          # Threshold calculator
search_motif.cpp     # Parallel motif search
README.md            # This file
```

---

## Configuration (`config.json`)

Edit `config.json` to suit your data and parameters:

```jsonc
{
  "paths": {
    "original_pwm_file": "path/to/input.meme",
    "sequence_file": "path/to/seq.fasta",
    "results_file": "results.txt"
  },
  "parameters": {
    "pwm2fm": 45, // Pseudocount multiplier
    "bayes": true, // Enable Bayesian smoothing
    "thr1": 95.0, // Full‐motif threshold percentile
    "thr2": 75.0, // Core‐region threshold percentile
    "simulation_iterations": 10000 // Monte Carlo iterations
  },
  "settings": {
    "core_length": 5, // Core window length for pre‐filtering
    "precision": 6 // Decimal places in output
  }
}
```

- **paths.original_pwm_file**: Input MEME motif file or 
- **paths.sequence_file**: Target sequences in FASTA format
- **paths.results_file**: Final hit list (TSV)
- **parameters.pwm2fm**: Pseudocount factor (or normalization constant)
- **parameters.bayes**: `true` to apply Bayesian smoothing, else simple normalization
- **parameters.thr1 / thr2**: Percentiles (0–100) for full and core scores
- **parameters.simulation_iterations**: More iterations → more stable thresholds, but slower
- **settings.core_length**: Short core region used to speed scanning
- **settings.precision**: Number of decimal places in PWM/threshold output

---

## Build

Make the scripts executable and compile:

```fish
chmod +x compile.sh run.sh
./compile.sh
```

Binaries produced:

- `meme2pwm`
- `get_thr`
- `search_motif`

---

## Run

Once compiled, simply run:

```fish
./run.sh
```

This will:

1. Convert MEME → PWM:
   ```
   ./meme2pwm <original_pwm> <TMP/pwm.txt> <bayes_flag> <pwm2fm>
   ```
2. Compute thresholds:
   ```
   ./get_thr <TMP/pwm.txt> <TMP/thr.txt> <thr1> <thr2> <iterations> <core_length> <precision>
   ```
3. Search motifs in parallel:
   ```
   ./search_motif <fasta> <TMP/pwm.txt> <TMP/thr.txt> <results_file> <core_length>
   ```

Intermediate files live in a temporary directory and are cleaned up automatically. The final TSV is saved to `results_file`.

---

## Interpreting Results

The output TSV (`results.txt` by default) has one hit per line, with columns:

| Column         | Description                                  |
| -------------- | -------------------------------------------- |
| **Sequence**   | FASTA header (e.g. `>chr1:100-500`)          |
| **MotifID**    | Motif identifier from the PWM file           |
| **Position**   | 0‐based start index in the strand            |
| **Strand**     | `positive` or `negative`                     |
| **Score**      | Log‐likelihood score (lower = better match)  |
| **MatchedSeq** | The matched DNA substring                    |


---

## License & Citation

Please cite this tool as “Goodness Of Alignment Test (GOAT), [Your Publication], 2025.”
