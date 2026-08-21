#!/usr/bin/env python3
"""Build a TElocal .locInd file from a TE GTF.

TElocal requires a pickled ``TElocal_Toolkit.TEindex.TEfeatures`` object
(the ``.locInd`` file).  This script builds one from a standard TE GTF
(plain or gzipped) using the same logic as the standalone
``TElocal_indexer`` binary, but avoids the external dependency by using
the ``TElocal_Toolkit`` library directly.

Usage:
  build_telocal_index.py --gtf te.gtf --out results/telocal/telocal.locInd
"""
import argparse
import os
import pickle
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read


def main():
    ap = argparse.ArgumentParser(
        description="Build a TElocal .locInd index from a TE GTF."
    )
    ap.add_argument("--gtf", required=True, help="TE GTF (plain or gzipped)")
    ap.add_argument("--out", required=True, help="Output .locInd path")
    args = ap.parse_args()

    # TEfeatures.build() requires a sorted, uncompressed GTF.
    # Sort by chr, then 1-based start position — matching the original
    # TElocal_indexer sort order.
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".gtf", delete=False
    ) as tmp:
        tmp_path = tmp.name
        with open_read(args.gtf) as fh:
            for line in fh:
                tmp.write(line)

    try:
        sorted_path = tmp_path + ".sorted"
        subprocess.check_call(
            ["sort", "-k1,1", "-k4,4n", tmp_path],
            stdout=open(sorted_path, "w"),
        )
        os.unlink(tmp_path)

        from TElocal_Toolkit.TEindex import TEfeatures

        te_idx = TEfeatures()
        te_idx.build(sorted_path)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        if os.path.exists(sorted_path):
            os.unlink(sorted_path)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "wb") as fh:
        pickle.dump(te_idx, fh)

    print(f"Built {args.out} with {te_idx.numInstances()} TE instances")


if __name__ == "__main__":
    main()
