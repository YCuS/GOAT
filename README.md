# GOAT

GOAT is a Linux-oriented PWM-based motif search pipeline for DNA FASTA files. It converts motif files into normalized PWMs, estimates fixed score thresholds, and scans sequences using a core-first strategy.

## Requirements

- Linux
- Bash
- C++17 compiler, such as `g++`

On Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y build-essential
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

```bash
bash compile.sh
bash run.sh examples/motifs.txt examples/sequences.fa results/motif_hits.tsv
```

Outputs:

- `results/thresholds.tsv`
- `results/motif_hits.tsv`

## Run

```bash
bash run.sh <motif_file> <sequence_file_or_dir> <results_file_or_dir> [options]
```

Required arguments:

- `motif_file`: MEME file or simple `Motif:` file.
- `sequence_file_or_dir`: FASTA file, or directory of `.fa`, `.fasta`, or `.fna` files.
- `results_file_or_dir`: output TSV file; output directory when sequence input is a directory.

Options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--threshold-file PATH` | next to results | Threshold output TSV. |
| `--core-length N` | `5` | Length of the selected core window. |
| `--full-pct PCT` | `95.0` | Full motif threshold percentile. |
| `--core-pct PCT` | `75.0` | Core region threshold percentile. |
| `--iterations N` | `10000` | Simulation count per motif. |
| `--precision N` | `6` | Threshold decimal precision. |
| `--bayes 0|1` | `0` | Bayesianize/smooth PWM probabilities in `meme2pwm`. |
| `--pseudocount N` | `45` | Pseudocount value passed to `meme2pwm`. |

`--bayes 1` enables Bayesian PWM smoothing during motif normalization. This reduces zero-probability issues by blending each PWM row with a small prior. `--pseudocount` controls how strongly the original motif row is weighted relative to that prior.

Example with options:

```bash
bash run.sh motifs.txt seq.fa results/hits.tsv \
  --threshold-file results/thresholds.tsv \
  --core-length 6 \
  --full-pct 95 \
  --core-pct 75 \
  --iterations 20000
```

## Batch Running

Batch mode is automatic when the sequence input is a directory:

```bash
bash run.sh motifs.txt fasta_dir results/batch_hits
```

The directory may contain files ending in `.fa`, `.fasta`, or `.fna`. GOAT writes one result file per FASTA input:

```text
results/batch_hits/<input_stem>_motif_hits.tsv
```

The threshold file is written once and reused for all FASTA files:

```text
results/batch_hits/thresholds.tsv
```

Batch mode also accepts the same threshold and PWM options:

```bash
bash run.sh motifs.txt fasta_dir results/batch_hits \
  --threshold-file results/batch_hits/thresholds.tsv \
  --core-length 6 \
  --bayes 1
```

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

Because scores are defined as negative log-likelihoods, lower values indicate better motif matches.

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
