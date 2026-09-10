#!/usr/bin/env python3
"""Cross-sample row totals for a feature x sample matrix, restricted to a
requested set of keys.

Used by chimera_candidates_explorer.R's rule (chimera_reads.smk) to join a
real cohort-total count onto each chimera candidate from
results/telocal/counts_matrix.tsv.gz (row key: TElocal locus key) and
results/chimera/assembly/counts_matrix.tsv.gz (row key: transcript_id),
without ever loading either matrix into R. Both can be genome/cohort-scale
(millions of loci for TElocal on a real annotation -- the same scale that
OOM'd telocal_counts once, see tecount_counts.py's sys.intern() fix) and
candidates_explorer.R already has its own prior OOM history from loading a
genome-wide reference file (te.bed) unfiltered. Streaming the matrix ONCE
here and emitting only the few thousand rows candidates.tsv.gz actually
references keeps R's input the same small size as its existing
gene_id_to_name.tsv.gz / genes.bed / te.bed joins.

Usage:
  chimera_candidates_matrix_totals.py --matrix counts_matrix.tsv.gz \\
      --keys candidate_keys.txt --out totals.tsv
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def load_keys(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def sum_matrix(matrix_path, wanted):
    """{key: total} for every row in `matrix_path` whose key is in `wanted`.
    Streams the file once; never holds more than `wanted`'s rows in memory."""
    totals = {}
    with open_read(matrix_path) as fh:
        fh.readline()  # header: key <TAB> sample1 <TAB> sample2 ...
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0]
            if key not in wanted:
                continue
            totals[key] = sum(float(v) for v in parts[1:])
    return totals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True,
                     help="feature x sample matrix (gzipped tsv)")
    ap.add_argument("--keys", required=True,
                     help="one row key per line -- only these are summed")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    wanted = load_keys(args.keys)
    totals = sum_matrix(args.matrix, wanted)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("key\ttotal\n")
        for key in sorted(totals):
            fh.write(f"{key}\t{totals[key]:.3f}\n")

    print(f"{len(totals)}/{len(wanted)} requested keys found in {args.matrix} "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
