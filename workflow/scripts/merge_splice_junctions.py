#!/usr/bin/env python3
"""Merge and filter multiple STAR SJ.out.tab files into one pooled junction
file for cohort-wide 2-pass mapping (STAR manual section 9.1).

Junctions are deduplicated by (chrom, start, end, strand), summing the
uniquely-mapping and multi-mapping read counts across all input files. Only
NOVEL junctions (STAR's own annotated flag == 0 in every input -- already-
annotated junctions are redundant since every pass already uses the shared,
annotation-guided genome index) with combined unique-read support >=
--min-unique-reads are kept, per the manual's recommendation to filter out
"likely false positives ... supported by a few reads" before reuse.

The output is written in the same 9-column SJ.out.tab format, which STAR's
--sjdbFileChrStartEnd accepts directly (manual section 9.1: "listing
SJ.out.tab files from all samples").

Usage:
  merge_splice_junctions.py --sj-tables a_SJ.out.tab b_SJ.out.tab ... \\
      --min-unique-reads 3 --out merged_SJ.out.tab
"""
import argparse


def main():
    ap = argparse.ArgumentParser(
        description="Merge/filter SJ.out.tab files for STAR cohort-wide 2-pass mapping."
    )
    ap.add_argument("--sj-tables", required=True, nargs="+")
    ap.add_argument("--min-unique-reads", type=int, default=3)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # key -> [motif, annotated, unique_reads, multi_reads, max_overhang]
    junctions = {}
    n_lines = 0
    for path in args.sj_tables:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                n_lines += 1
                chrom, start, end, strand = cols[0], cols[1], cols[2], cols[3]
                motif, annotated = cols[4], cols[5]
                unique_reads, multi_reads, overhang = (
                    int(cols[6]), int(cols[7]), int(cols[8]),
                )
                key = (chrom, int(start), int(end), strand)
                if key not in junctions:
                    junctions[key] = [motif, annotated, 0, 0, 0]
                j = junctions[key]
                j[2] += unique_reads
                j[3] += multi_reads
                j[4] = max(j[4], overhang)

    kept = []
    for (chrom, start, end, strand), (motif, annotated, uniq, multi, overhang) in junctions.items():
        if annotated != "0":
            continue  # already-annotated junction: redundant, every pass already has it
        if uniq < args.min_unique_reads:
            continue
        kept.append((chrom, start, end, strand, motif, annotated, uniq, multi, overhang))

    kept.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(args.out, "w") as fh:
        for row in kept:
            fh.write("\t".join(str(v) for v in row) + "\n")

    print(
        f"{n_lines} junction records read across {len(args.sj_tables)} sample(s), "
        f"{len(junctions)} unique junctions, {len(kept)} kept "
        f"(novel, >= {args.min_unique_reads} combined unique reads)"
    )


if __name__ == "__main__":
    main()
