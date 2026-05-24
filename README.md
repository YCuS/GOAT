# GOAT

GOAT is a Linux-oriented PWM motif search pipeline for DNA FASTA files. It converts motif files to normalized PWMs, estimates fixed score thresholds, and searches sequences with a core-first scan.

## Requirements

- Linux
- Bash
- `jq`
- C++17 compiler, such as `g++`

On Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y build-essential jq
```

## Build

```bash
bash compile.sh
```

This creates:

- `meme2pwm`
- `get_thr`
- `search_motif`

## Quick Start

The default `config.json` uses the example files in `examples/`.

```bash
bash compile.sh
bash run.sh
```

Outputs:

- `results/thresholds.tsv`
- `results/motif_hits.tsv`

## Configuration

Edit `config.json` for your data:

```json
{
  "parameters": {
    "pwm2fm": 45,
    "bayes": false,
    "thr1": 95.0,
    "thr2": 75.0,
    "simulation_iterations": 10000
  },
  "paths": {
    "original_pwm_file": "examples/motifs.txt",
    "sequence_file": "examples/sequences.fa",
    "results_file": "results/motif_hits.tsv",
    "threshold_file": "results/thresholds.tsv"
  },
  "settings": {
    "save_threshold": true,
    "core_length": 5,
    "precision": 6
  }
}
```

Key fields:

- `paths.original_pwm_file`: motif input file, either MEME format or simple `Motif:` format.
- `paths.sequence_file`: FASTA file, or a directory of `.fa`, `.fasta`, or `.fna` files.
- `paths.results_file`: output TSV file, or output directory in batch mode.
- `paths.threshold_file`: threshold output TSV.
- `parameters.thr1`: full motif threshold percentile.
- `parameters.thr2`: core region threshold percentile.
- `parameters.simulation_iterations`: simulation count per motif.
- `settings.core_length`: length of the selected core window.

## Motif Input

Simple motif format:

```text
Motif:example_motif
0.90 0.03 0.04 0.03
0.03 0.90 0.04 0.03
0.03 0.04 0.90 0.03
0.03 0.04 0.03 0.90
```

Each row is ordered as `A C G T`. Rows are normalized automatically. MEME files are also supported when they contain `letter-probability matrix:`.

## Search Method

For each motif, `get_thr` selects the lowest-entropy contiguous core window. During search, `search_motif` first scores the core region at each candidate position. Only windows that pass the core threshold are scored against the full motif. Sequence scanning is parallelized across available CPU threads.

Scores are negative log-likelihoods, so lower values indicate better matches.

## Output

Threshold file:

```text
# motif_id	full_threshold	core_threshold	core_start	motif_length
```

Motif hit file:

```text
sequence_id	motif_id	strand	start	end	full_score	core_score	matched_sequence
```

Coordinates are 1-based and inclusive on the input FASTA sequence.

## Manual Usage

```bash
./meme2pwm examples/motifs.txt /tmp/pwm.tsv 0 45
./get_thr /tmp/pwm.tsv results/thresholds.tsv 95.0 75.0 10000 5 6
./search_motif examples/sequences.fa /tmp/pwm.tsv results/thresholds.tsv results/motif_hits.tsv 5
```

Batch mode:

```bash
./search_motif path/to/fasta_dir /tmp/pwm.tsv results/thresholds.tsv results/batch_hits 5 --batch
```
