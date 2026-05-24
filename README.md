# GOAT

GOAT is a lightweight PWM-only motif search pipeline for DNA FASTA files. It keeps the workflow intentionally direct:

1. Convert motifs into normalized position weight matrices.
2. Simulate fixed score thresholds from each PWM.
3. Search sequences with a core-first filter, then score the full motif only for windows that pass the core threshold.

This repository does not use GC-content adjustment, motif clustering, relaxed/strict cluster rules, or gap alignment.

## Features

- Supports MEME motif files and a simple `Motif:` PWM-like input format.
- Normalizes motif rows automatically as `A C G T`.
- Selects the motif core as the lowest-entropy contiguous window.
- Uses Monte Carlo simulation to calculate fixed thresholds per motif.
- Searches both `+` and `-` strands.
- Applies a fast core gate before full motif scoring.
- Uses multi-threaded sequence scanning for faster motif search.
- Writes tab-separated outputs with stable headers.
- Supports single FASTA input or batch mode over a directory of FASTA files.

## Requirements

- Bash
- `jq`
- A C++17 compiler, such as `g++` or `clang++`

On macOS with Homebrew:

```bash
brew install jq
```

## Build

```bash
bash compile.sh
```

This creates:

- `meme2pwm`: converts MEME or simple motif input into normalized PWM format
- `get_thr`: simulates fixed full-motif and core-region thresholds
- `search_motif`: scans FASTA sequences and writes motif hits

## Quick Start

The default `config.json` uses the files under `examples/`.

```bash
bash compile.sh
bash run.sh
```

Expected outputs:

- `results/thresholds.tsv`
- `results/motif_hits.tsv`

## Configuration

Edit `config.json` before running on your own data:

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

### Parameters

| Field | Meaning |
| --- | --- |
| `parameters.pwm2fm` | Pseudocount-related value passed to `meme2pwm` when Bayesian smoothing is enabled. |
| `parameters.bayes` | Whether to apply Bayesian smoothing during PWM normalization. |
| `parameters.thr1` | Full-motif score percentile used as the full threshold. |
| `parameters.thr2` | Core-region score percentile used as the core threshold. |
| `parameters.simulation_iterations` | Number of simulated motif sequences per motif. |
| `settings.core_length` | Length of the contiguous core window. |
| `settings.precision` | Decimal precision for threshold output. |

### Paths

| Field | Meaning |
| --- | --- |
| `paths.original_pwm_file` | Motif input in MEME format or simple `Motif:` format. |
| `paths.sequence_file` | FASTA file, or a directory containing `.fa`, `.fasta`, or `.fna` files. |
| `paths.results_file` | Output TSV file for single input; output directory in batch mode. |
| `paths.threshold_file` | Output threshold TSV file. |

## Motif Input

### Simple Format

Each motif starts with `Motif:<id>`. Every following row is one PWM position in `A C G T` order.

```text
Motif:example_motif
0.90 0.03 0.04 0.03
0.03 0.90 0.04 0.03
0.03 0.04 0.90 0.03
0.03 0.04 0.03 0.90
```

Rows are normalized automatically, so counts or probabilities are both acceptable as long as the columns are ordered as `A C G T`.

### MEME Format

Files containing `letter-probability matrix:` are detected as MEME format and parsed by `meme2pwm`.

## Algorithm

### 1. PWM Normalization

`meme2pwm` reads motifs and writes a normalized internal PWM file:

```text
Motif:<motif_id>
A_probability  C_probability  G_probability  T_probability
...
```

### 2. Core Selection

For each motif, `get_thr` selects a contiguous core window of length `settings.core_length`.

The chosen core is the window with the smallest entropy sum across its columns. In practice, this favors the most informative and least ambiguous region of the PWM.

### 3. Threshold Simulation

For each motif:

- Simulated sequences are sampled from the motif PWM.
- Each simulated sequence is scored against the full PWM.
- The selected core region is scored separately.
- `thr1` and `thr2` percentiles become fixed full and core thresholds.

Scores are negative log-likelihoods. Lower values indicate better matches.

### 4. Core-First Motif Search

`search_motif` scans each candidate window in this order:

1. Score only the core region at the expected core offset.
2. Reject the window immediately if the core score is above the core threshold.
3. Score the full motif only after the core passes.
4. Report the hit only if the full score also passes the full threshold.

This keeps the original PWM-only logic while avoiding unnecessary full-motif scoring for weak windows.

### 5. Parallel Search

Sequence scanning is parallelized across available CPU threads. Each worker scans a chunk of sequences independently, then the results are merged and sorted for stable output.

## Output

### Thresholds

`thresholds.tsv`:

```text
# motif_id	full_threshold	core_threshold	core_start	motif_length
```

| Column | Meaning |
| --- | --- |
| `motif_id` | Motif identifier. |
| `full_threshold` | Maximum accepted full-motif negative log-likelihood score. |
| `core_threshold` | Maximum accepted core-region negative log-likelihood score. |
| `core_start` | 0-based core start position inside the motif. |
| `motif_length` | Motif length in bases. |

### Motif Hits

`motif_hits.tsv`:

```text
sequence_id	motif_id	strand	start	end	full_score	core_score	matched_sequence
```

| Column | Meaning |
| --- | --- |
| `sequence_id` | FASTA record identifier. |
| `motif_id` | Matched motif identifier. |
| `strand` | `+` or `-`. |
| `start` | 1-based inclusive start coordinate on the input sequence. |
| `end` | 1-based inclusive end coordinate on the input sequence. |
| `full_score` | Full-motif negative log-likelihood score. |
| `core_score` | Core-region negative log-likelihood score. |
| `matched_sequence` | Matched sequence in motif orientation. |

## Manual Commands

The pipeline can also be run step by step:

```bash
./meme2pwm examples/motifs.txt /tmp/pwm.tsv 0 45
./get_thr /tmp/pwm.tsv results/thresholds.tsv 95.0 75.0 10000 5 6
./search_motif examples/sequences.fa /tmp/pwm.tsv results/thresholds.tsv results/motif_hits.tsv 5
```

For batch input:

```bash
./search_motif path/to/fasta_dir /tmp/pwm.tsv results/thresholds.tsv results/batch_hits 5 --batch
```

## Notes

- Coordinates are reported on the original FASTA sequence.
- `matched_sequence` for `-` strand hits is shown in motif orientation.
- The current implementation is intentionally PWM-only and fixed-threshold based.
- Generated binaries and `results/` are ignored by git.
