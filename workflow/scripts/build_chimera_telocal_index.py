#!/usr/bin/env python3
"""Build the merged, all-samples TElocal interval index used by
chimera_telocal_annotate.py.

This runs ONCE per pipeline run (chimera_telocal_index rule) instead of once
per sample: chimera_telocal_annotate previously re-parsed and re-summed every
sample's TElocal cntTable from scratch in every one of its N per-sample jobs,
an Nx redundant rebuild of the identical shared index. Building it once here
and having every annotate job just deserialize the (compact, columnar --
see chimera_telocal_index.py) result removes that redundant work and lets
the one expensive build be sized on its own, instead of forcing every
lightweight per-sample annotate job to pay for it.

Usage:
  build_chimera_telocal_index.py \\
      --telocal-tables {sample1}.cntTable.gz {sample2}.cntTable.gz ... \\
      --out telocal_index.pkl.gz
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chimera_telocal_index import TelocalIndex
from gz_io import open_read


def main():
    ap = argparse.ArgumentParser(
        description="Build the merged TElocal interval index for the chimera screen."
    )
    ap.add_argument("--telocal-tables", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    index = TelocalIndex.build(args.telocal_tables, open_read)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    index.save(args.out)

    n_loci = sum(len(cols["key"]) for cols in index.chroms.values())
    print(f"{n_loci} TE loci across {len(index.chroms)} chromosomes from "
          f"{len(args.telocal_tables)} sample tables -> {args.out}")


if __name__ == "__main__":
    main()
