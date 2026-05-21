#!/usr/bin/env python3
"""
transfac2pwm.py

Convert a JASPAR / TRANSFAC style combined matrix file into the simplified PWM
format used by this project (the same output style produced by meme2pwm):

Output format per motif:
    Motif:<motif_id>
    pA\tpC\tpG\tpT   (one line per position, probabilities, tab separated)
    <blank line>

Supported input patterns (auto-detected per motif):
1) JASPAR PFM style (four lines per motif, starting with A / C / G / T followed by counts columns)
   Example:
       >MA0001.1 AGL3
       A  1  5  2  0
       C  0  2  3  4
       G  3  1  0  2
       T  4  0  3  1

2) TRANSFAC-like positional rows (position index + four counts) *[optional]*
   Example:
       ID MA0001.1
       01  1 5 2 0
       02  0 2 3 4
       03  3 1 0 2
       04  4 0 3 1
       //

3) Simple already pivoted format (A C G T columns separated by whitespace) preceded by a header line
   starting with 'Motif:' (these will simply be passed through after normalization if requested).

Probabilities are computed with optional pseudocount smoothing:
    prob = (count + pseudocount) / (row_sum + 4 * pseudocount)

If the input presents counts per base (A-lines, C-lines, etc.), counts are pivoted into per-position rows.

Usage:
    ./transfac2pwm.py -i input.txt -o output_pwm.txt [--pseudocount 0.5] [--min-len 2]

Arguments:
    -i / --input         Input matrix file
    -o / --output        Output PWM file
    --pseudocount FLOAT  (default: 0.0)
    --min-len INT        Skip motifs shorter than this length (default: 1)
    --uppercase-id       Force motif IDs to uppercase
    --keep-empty         Do not drop motifs that would otherwise be skipped (diagnostic)

Exit codes:
    0 success
    1 command line / file error
    2 no motifs written

Author: Auto-generated helper script
"""
from __future__ import annotations
import sys
import argparse
import re
from typing import List, Dict, Optional

BASES = ["A", "C", "G", "T"]
BASE_SET = set(BASES)

motif_header_re = re.compile(r"^>(\S+)")  # >MotifID ...
transfac_id_re = re.compile(r"^ID\s+(\S+)")  # ID motif_id
transfac_row_re = re.compile(
    r"^(\d+)[\t ]+([\d.eE+-]+)[\t ]+([\d.eE+-]+)[\t ]+([\d.eE+-]+)[\t ]+([\d.eE+-]+)\s*$"
)
simple_pwm_row_re = re.compile(r"^[\d.eE+-]+(\s+[\d.eE+-]+){3}\s*$")

# Lines like: A  1 2 3 4 5 ...
base_line_re = re.compile(r"^([ACGTacgt])(?:\s+|\t)([\d.eE+\-\s\t]+)$")


def parse_args(argv: List[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Convert JASPAR/TRANSFAC matrices to project PWM format"
    )
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input matrix file (JASPAR/TRANSFAC/combined)",
    )
    ap.add_argument("-o", "--output", required=True, help="Output PWM file")
    ap.add_argument(
        "--pseudocount",
        type=float,
        default=0.0,
        help="Pseudocount added to each cell (default 0)",
    )
    ap.add_argument(
        "--min-len",
        type=int,
        default=1,
        help="Minimum motif length to keep (default 1)",
    )
    ap.add_argument(
        "--uppercase-id", action="store_true", help="Force motif IDs to uppercase"
    )
    ap.add_argument(
        "--keep-empty",
        action="store_true",
        help="Write motifs even if they would be skipped (diagnostic)",
    )
    return ap.parse_args(argv)


class MotifBuilder:
    def __init__(self):
        self.motif_id: Optional[str] = None
        self.counts_rows: List[List[float]] = (
            []
        )  # List of [A,C,G,T] counts per position (pos-based)
        self.base_counts: Dict[str, List[float]] = {}  # For A/C/G/T lines form
        self.mode: Optional[str] = None  # 'pfm_base', 'pos_rows', 'simple', or None

    def start_new(self, motif_id: str):
        self.finalize()  # flush existing
        self.motif_id = motif_id
        self.counts_rows = []
        self.base_counts = {}
        self.mode = None

    def add_base_line(self, base: str, numbers: List[float]):
        base = base.upper()
        if self.mode is None:
            self.mode = "pfm_base"
        elif self.mode != "pfm_base":
            # Mixed mode, ignore line
            return
        self.base_counts[base] = numbers

    def add_position_row(self, row: List[float]):
        if self.mode is None:
            self.mode = "pos_rows"
        if self.mode == "pos_rows":
            self.counts_rows.append(row)

    def add_simple_row(self, row: List[float]):
        if self.mode is None:
            self.mode = "simple"
        if self.mode == "simple":
            self.counts_rows.append(row)

    def finalize(self) -> Optional[Dict]:
        if not self.motif_id:
            return None
        # Convert base counts (A,C,G,T lines) into position rows
        if self.mode == "pfm_base":
            lengths = {b: len(v) for b, v in self.base_counts.items()}
            if len(lengths) == 4 and len(set(lengths.values())) == 1:
                L = next(iter(lengths.values()))
                self.counts_rows = []
                for i in range(L):
                    self.counts_rows.append(
                        [
                            self.base_counts.get("A", [0] * L)[i],
                            self.base_counts.get("C", [0] * L)[i],
                            self.base_counts.get("G", [0] * L)[i],
                            self.base_counts.get("T", [0] * L)[i],
                        ]
                    )
            else:
                # inconsistent lengths, drop
                data = self._dump_and_reset()
                return None
        data = None
        if self.counts_rows:
            data = {
                "id": self.motif_id,
                "counts": self.counts_rows,  # list of [A,C,G,T]
            }
        # Reset for next motif
        self.motif_id = None
        self.counts_rows = []
        self.base_counts = {}
        self.mode = None
        return data

    def _dump_and_reset(self):
        self.motif_id = None
        self.counts_rows = []
        self.base_counts = {}
        self.mode = None


def iterate_motifs(path: str):
    builder = MotifBuilder()
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()  # keep simple trimming
            if not line:
                # Blank line: finalize motif (for PFM grouped style)
                m = builder.finalize()
                if m:
                    yield m
                continue

            if line.startswith("//"):
                m = builder.finalize()
                if m:
                    yield m
                continue

            # Header patterns
            m = motif_header_re.match(line)
            if m:
                builder.start_new(m.group(1))
                continue
            m = transfac_id_re.match(line)
            if m:
                builder.start_new(m.group(1))
                continue

            # Base lines (A/C/G/T counts across positions)
            mb = base_line_re.match(line)
            if mb:
                base = mb.group(1)
                nums = [float(x) for x in mb.group(2).split() if x]
                builder.add_base_line(base, nums)
                continue

            # TRANSFAC positional row
            mt = transfac_row_re.match(line)
            if mt:
                row = [float(mt.group(i)) for i in range(2, 6)]
                builder.add_position_row(row)
                continue

            # Simple row (A C G T already) if we have a motif active
            if builder.motif_id and simple_pwm_row_re.match(line):
                row = [float(x) for x in line.split()]
                if len(row) == 4:
                    builder.add_simple_row(row)
                continue
        # EOF flush
        m = builder.finalize()
        if m:
            yield m


def counts_to_probabilities(
    count_rows: List[List[float]], pseudocount: float
) -> List[List[float]]:
    probs: List[List[float]] = []
    for row in count_rows:
        if len(row) != 4:
            continue
        row_sum = sum(row) + 4.0 * pseudocount
        if row_sum <= 0:
            probs.append([0.25, 0.25, 0.25, 0.25])
            continue
        probs.append([(c + pseudocount) / row_sum for c in row])
    return probs


def write_pwm(
    motifs: List[Dict],
    out_path: str,
    pseudocount: float,
    min_len: int,
    uppercase: bool,
    keep_empty: bool,
) -> int:
    written = 0
    with open(out_path, "w") as out:
        for m in motifs:
            motif_id = m["id"]
            if uppercase:
                motif_id = motif_id.upper()
            probs = counts_to_probabilities(m["counts"], pseudocount)
            if not probs or len(probs) < min_len:
                if not keep_empty:
                    continue
            out.write(f"Motif:{motif_id}\n")
            for row in probs:
                out.write("\t".join(f"{p:.6f}" for p in row) + "\n")
            out.write("\n")
            written += 1
    return written


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    try:
        motifs = list(iterate_motifs(args.input))
        if not motifs:
            print("No motifs parsed from input", file=sys.stderr)
        written = write_pwm(
            motifs,
            args.output,
            args.pseudocount,
            args.min_len,
            args.uppercase_id,
            args.keep_empty,
        )
        if written == 0:
            print("No motifs written (after filtering)", file=sys.stderr)
            return 2
        print(f"Converted {written} motifs -> {args.output}")
        return 0
    except FileNotFoundError:
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
