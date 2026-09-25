#!/usr/bin/env python3
"""Merge the per-sample SJ.out.tab-derived junction tables
(classify_chimera_splice_junctions.py) into the all-events catalog and the event x sample
counts matrix -- the SJ-based screen's counterpart of chimera_chimeric_reads_counts.py,
same shape, same conventions.

Outputs:

  all_events.tsv   one row per unique event_id across all samples, with the
                   annotation columns taken from the first sample that saw it
                   (event_id is breakpoint-deterministic, so annotations
                   agree across samples) plus n_samples / total_reads.
  counts_matrix.tsv  event_id x sample matrix of STAR's own unique_reads
                   count for that junction (0 where a sample never had this
                   junction). Written when --out-counts is given.
  cpm_matrix.tsv   the same matrix, each column divided by that sample's
                   column total x 1e6 (CPM, not TPM -- a splice junction has
                   no meaningful "length" to normalize by, same rationale as
                   chimera_chimeric_reads_counts.py). Written when --out-cpm is given.
  te-gene-junctions.tsv  the all_events catalog filtered to gene<->TE events
                   (direction gene_to_te / te_to_gene), written when
                   --out-te-events is given.

Nothing is filtered here beyond what classify_chimera_splice_junctions.py itself already
applied (--min-unique-reads / --require-canonical) -- annotate-first, same as
the other two screens.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

ANNOTATION_COLUMNS = [
    "event_id", "chrom", "intron_start", "intron_end", "strand", "motif",
    "canonical", "annotated", "multi_reads", "overhang",
    "donor_hits", "acceptor_hits", "direction", "direction_ambiguous",
    "gene_id", "gene_strand", "te_id", "te_subfamily", "te_family",
    "te_class", "chimera_type", "antisense_flag", "library_strand",
    "transcript_strand", "gene_strand_match",
]


def load(path):
    rows = []
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--sample-names", required=True, nargs="+")
    ap.add_argument("--out-events", required=True)
    ap.add_argument("--out-counts", required=False)
    ap.add_argument("--out-cpm", required=False)
    ap.add_argument(
        "--out-te-events", required=False,
        help="Optional output: the all-events catalog filtered to gene<->TE "
        "events (direction gene_to_te / te_to_gene).",
    )
    args = ap.parse_args()

    if len(args.tables) != len(args.sample_names):
        sys.exit("error: --tables and --sample-names must have equal length")

    events = {}
    for path, sample in zip(args.tables, args.sample_names):
        for row in load(path):
            eid = row["event_id"]
            ev = events.setdefault(eid, {"sample": sample, "counts": {}})
            ev["counts"][sample] = int(row["unique_reads"])
            if ev["sample"] == sample:
                for col in ANNOTATION_COLUMNS:
                    ev[col] = row.get(col, ".")
    order = sorted(events, key=lambda e: (events[e]["sample"], e))

    with open_write(args.out_events) as fh:
        header = ANNOTATION_COLUMNS + ["n_samples", "total_reads"]
        fh.write("\t".join(header) + "\n")
        for eid in order:
            ev = events[eid]
            row = [ev.get(c, ".") for c in ANNOTATION_COLUMNS]
            row += [len(ev["counts"]), sum(ev["counts"].values())]
            fh.write("\t".join(str(x) for x in row) + "\n")

    if args.out_counts:
        with open_write(args.out_counts) as fh:
            fh.write("event_id\t" + "\t".join(args.sample_names) + "\n")
            for eid in order:
                counts = events[eid]["counts"]
                vals = [counts.get(s, 0) for s in args.sample_names]
                fh.write(eid + "\t" + "\t".join(str(c) for c in vals) + "\n")

    if args.out_cpm:
        totals = dict.fromkeys(args.sample_names, 0)
        for eid in order:
            counts = events[eid]["counts"]
            for s in args.sample_names:
                totals[s] += counts.get(s, 0)
        with open_write(args.out_cpm) as fh:
            fh.write("event_id\t" + "\t".join(args.sample_names) + "\n")
            for eid in order:
                counts = events[eid]["counts"]
                cpm = [
                    0 if totals[s] == 0 else counts.get(s, 0) / totals[s] * 1e6
                    for s in args.sample_names
                ]
                fh.write(eid + "\t" + "\t".join(f"{v:.3f}" for v in cpm) + "\n")

    if args.out_te_events:
        with open_write(args.out_te_events) as fh:
            fh.write("\t".join(ANNOTATION_COLUMNS + ["n_samples", "total_reads"]) + "\n")
            for eid in order:
                ev = events[eid]
                if ev.get("direction") not in ("gene_to_te", "te_to_gene"):
                    continue
                row = [ev.get(c, ".") for c in ANNOTATION_COLUMNS]
                row += [len(ev["counts"]), sum(ev["counts"].values())]
                fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{len(order)} unique junctions across {len(args.sample_names)} samples "
          f"-> {args.out_events}")


if __name__ == "__main__":
    main()
