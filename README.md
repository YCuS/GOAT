# Goodness Of Alignment Test (GOAT)

This is a tool for searching motifs in DNA/RNA sequences using Position Weight Matrices (PWMs) and Monte Carlo simulation, optimized for Linux systems.

---

## Prerequisites

- Linux
- A C++17 compiler (e.g. `g++`)
- `jq` (for JSON parsing)
- Bash shell

---

## Repository Structure

```bash
compile.sh           # Build all binaries
run.sh               # Orchestrate config.json workflow
config.json          # User parameters and paths
meme2pwm.cpp         # MEME→PWM converter
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
    "sequence_file": "path/to/input.fasta or directory", // FASTA file or directory of FASTA files
    "results_file": "results.txt or directory" // Output file (single) or directory (batch)
  },
  "parameters": {
    "pwm2fm": 45, // Pseudocount multiplier
    "bayes": true, // Enable Bayesian smoothing
    "thr1": 95.0, // Full-motif threshold percentile
    "thr2": 75.0, // Core-region threshold percentile
    "simulation_iterations": 10000 // Monte Carlo iterations
  },
  "settings": {
    "core_length": 5, // Core window length for pre-filtering
    "precision": 6 // Decimal places in output
  }
}
```

- **paths.original_pwm_file**: Input MEME motif file
- **paths.sequence_file**: FASTA file or directory of FASTA files for batch mode
- **paths.results_file**: Output TSV file (single) or directory (batch mode)

---

## Build

Make the scripts executable and compile:

```bash
chmod +x compile.sh run.sh
./compile.sh
```

Binaries produced:

- `meme2pwm`
- `get_thr`
- `search_motif`

---

## Run

After building, configure `config.json` and execute:

```bash
chmod +x run.sh
./run.sh
```

`run.sh` auto-detects single-file vs. directory input:

- **Single-file mode**: `sequence_file` → one output file `results_file`
- **Batch mode**: `sequence_file` is a directory → `results_file` treated as output directory; generates one `<stem>_results.txt` per FASTA

Intermediate files are stored in a temporary directory and cleaned up automatically.

---

## Interpreting Results

The output TSV (default `results.txt`) has columns:

| Column     | Description                                 |
| ---------- | ------------------------------------------- |
| Sequence   | FASTA header (e.g. `>chr1:100-500`)         |
| MotifID    | Motif identifier                            |
| Position   | 0-based start index in the strand           |
| Strand     | `positive` or `negative`                    |
| Score      | Log-likelihood score (lower = better match) |
| MatchedSeq | The matched DNA substring                   |

---

## License & Citation

Please cite this tool as “Goodness Of Alignment Test (GOAT), [Your Publication], 2025.”
