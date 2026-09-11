#!/usr/bin/env python3
"""Sum results/chimera/assembly/counts_matrix.tsv.gz (transcript_id x sample,
StringTie-derived raw count ESTIMATES, see quantify_chimera_assembly.py) down
to one row per (matched_gene_id, te_id, chimera_type) -- the grain DESeq2/
edgeR should actually test at, not the raw MSTRG transcript_id.

Why this grain and not transcript_id or (gene, te) alone:

  * transcript_id is too fragmented for a DE test. StringTie routinely
    assembles one real transcriptional event as several near-identical
    isoforms differing by a few bp at an imprecise exon boundary (see
    classify_chimera_assembly.py's breakpoint-tolerance discussion) --
    splitting the read support for one event across N noisy rows.
  * Summing every isoform of a (gene_id, te_id) pair together, ignoring
    chimera_type, is wrong: a single candidate pair can carry isoforms of
    DIFFERENT chimera_type (te_initiated vs te_terminated vs te_exonized),
    which are mechanistically distinct regulatory events -- pooling them
    can cancel out or mask a real signal specific to just one mechanism.
  * (gene_id, te_id, chimera_type) only pools the sample dimension (never
    collapsed) and only pools isoforms that already share the same
    structural class -- the legitimate unit.

Rows whose transcripts.tsv.gz candidate has matched_gene_id=="." and/or
te_id=="." (te_initiated_intergenic has no matched gene; unspliced_te_only
has no splice evidence) are excluded, same filter idiom chimera_evidence.py
uses for its own (gene, te) pair aggregation.

Output (gene_te_chimera_id x sample, DESeq2/edgeR-ready):
  gene_te_chimera_id   f"{matched_gene_id}:{te_id}:{chimera_type}" -- a
                       single joined key column (not three separate
                       columns) so read.delim(file, row.names=1) in R
                       hands DESeqDataSetFromMatrix a pure numeric matrix
                       with no further munging. Split on ":" to recover
                       the three parts individually.
  {sample columns}     summed raw counts (integers), one column per sample,
                       in counts_matrix.tsv.gz's own header order.

Usage:
  aggregate_chimera_assembly_counts.py \\
      --transcripts results/chimera/assembly/transcripts.tsv.gz \\
      --counts results/chimera/assembly/counts_matrix.tsv.gz \\
      --out gene_te_chimera_counts_matrix.tsv.gz
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def load_transcript_groups(path):
    """{transcript_id: (gene, te, chimera_type) or None} -- None for rows
    that are not real gene-TE candidates (matched_gene_id or te_id == "."),
    same exclusion chimera_evidence.py applies. Every transcript_id is
    still a dict key (even when its value is None) so the counts-matrix
    pass below can tell "known, filtered out" apart from "unknown id",
    which would mean the two input files are out of sync."""
    groups = {}
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            tid = cols[idx["transcript_id"]]
            gene = cols[idx["matched_gene_id"]]
            te = cols[idx["te_id"]]
            ctype = cols[idx["chimera_type"]]
            if gene in (".", "") or te in (".", "") or ctype in (".", ""):
                groups[tid] = None
            else:
                groups[tid] = (gene, te, ctype)
    return groups


def aggregate(counts_path, groups):
    """{(gene, te, chimera_type): [per-sample summed int]}, plus the sample
    list (from counts_matrix.tsv.gz's own header) and include/exclude
    counts for the summary line."""
    sums = {}
    n_included = n_excluded = 0
    with open_read(counts_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            tid = parts[0]
            if tid not in groups:
                sys.exit(f"error: transcript_id {tid!r} in {counts_path} not "
                         f"found in --transcripts -- inputs out of sync")
            key = groups[tid]
            if key is None:
                n_excluded += 1
                continue
            n_included += 1
            values = [int(v) for v in parts[1:]]
            row = sums.setdefault(key, [0] * len(samples))
            for i, v in enumerate(values):
                row[i] += v
    return samples, sums, n_included, n_excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcripts", required=True,
                     help="results/chimera/assembly/transcripts.tsv.gz")
    ap.add_argument("--counts", required=True,
                     help="results/chimera/assembly/counts_matrix.tsv.gz")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    groups = load_transcript_groups(args.transcripts)
    samples, sums, n_included, n_excluded = aggregate(args.counts, groups)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open_write(args.out) as fh:
        fh.write("gene_te_chimera_id\t" + "\t".join(samples) + "\n")
        for gene, te, ctype in sorted(sums):
            row = sums[(gene, te, ctype)]
            fh.write(f"{gene}:{te}:{ctype}\t" + "\t".join(str(v) for v in row) + "\n")

    print(f"{len(sums)} (gene,te,chimera_type) groups from {n_included} "
          f"candidate transcripts ({n_excluded} excluded: no matched gene "
          f"and/or TE) -> {args.out}")


if __name__ == "__main__":
    main()
