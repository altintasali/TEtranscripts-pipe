#!/usr/bin/env python3
"""Cross-check assembly-based chimera candidates against the junction
screen's read-level gene<->TE calls (results/chimera_junction/te-gene-chimeras.tsv.gz).

The two screens use independent evidence (StringTie transcript structure vs.
STAR chimeric-junction reads), so a (gene_id, te_id) pair called by BOTH is
materially higher confidence than either alone -- this just adds a boolean
flag + the supporting junction-read count, it doesn't drop anything from
either table.
"""
import argparse
import gzip
import os


def read_table(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return header, rows


def open_write(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True,
                     help="results/chimera_assembly/candidates.tsv.gz")
    ap.add_argument("--junction", required=True,
                     help="results/chimera_junction/te-gene-chimeras.tsv.gz")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cand_header, cand_rows = read_table(args.candidates)
    junc_header, junc_rows = read_table(args.junction)

    junc_gene_idx = junc_header.index("gene_id")
    junc_te_idx = junc_header.index("te_id")
    junc_reads_idx = junc_header.index("reads")

    # (gene_id, te_id) -> total supporting reads across all junction events
    junc_support = {}
    for row in junc_rows:
        key = (row[junc_gene_idx], row[junc_te_idx])
        junc_support[key] = junc_support.get(key, 0) + int(row[junc_reads_idx])

    gene_idx = cand_header.index("matched_gene_id")
    te_idx = cand_header.index("te_id")

    out_header = cand_header + ["confirmed_by_junction_screen", "junction_supporting_reads"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("\t".join(out_header) + "\n")
        n_confirmed = 0
        for row in cand_rows:
            key = (row[gene_idx], row[te_idx])
            reads = junc_support.get(key)
            confirmed = "yes" if reads else "no"
            n_confirmed += confirmed == "yes"
            fh.write("\t".join(row + [confirmed, str(reads or 0)]) + "\n")

    print(f"{n_confirmed}/{len(cand_rows)} assembly-based candidates also have "
          f"chimeric-junction read support")


if __name__ == "__main__":
    main()
