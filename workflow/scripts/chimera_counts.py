#!/usr/bin/env python3
"""Merge the per-sample junction tables (classify_chimera_junctions.py) into the
chimera all-events catalog and the event x sample counts matrix.

Outputs:

  all_events.tsv   one row per unique event_id across all samples, with the
                   annotation columns taken from the first sample that saw it
                   (annotations are breakpoint-deterministic, so they agree
                   across samples) plus n_samples / total_reads.
  counts_matrix.tsv  event_id x sample read-count matrix (0 where a sample has
                   no reads supporting the event). Written only when the
                   chimera.outputs.write_counts_matrix config is true.
  te-gene-chimeras.tsv  the all_events catalog filtered to gene<->TE events
                   (direction gene_to_te / te_to_gene), written when
                   --out-te-events is given -- the TE chimeras as their own
                   table.

Nothing is filtered from the main outputs: the full union of annotated events
is shipped (the QC filters in the chimera.qc config section apply only to the
PCA/clustering view in sample_qc.smk).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

ANNOTATION_COLUMNS = [
    "event_id", "donor_chrom", "donor_breakpoint", "donor_strand",
    "acceptor_chrom", "acceptor_breakpoint", "acceptor_strand",
    "junction_type", "canonical",
    "repeat_flag", "direction", "direction_ambiguous",
    "gene_id", "gene_strand", "te_id", "te_subfamily",
    "te_family", "te_class", "chimera_type", "antisense_flag",
    "library_strand", "transcript_strand", "gene_strand_match",
    "telocal_count", "telocal_locus", "telocal_active",
    "te_refined_by_telocal",
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
            ev["counts"][sample] = int(row["reads"])
            if ev["sample"] == sample:
                # first-seen annotations (event_id is breakpoint-deterministic).
                # .get(): the telocal_* columns only exist when telocal is
                # enabled (chimera_counts_input() then feeds the
                # _with-telocal tables); with telocal off the plain junction
                # tables legitimately lack them, and strict indexing here
                # crashed the whole merge.
                for col in ANNOTATION_COLUMNS:
                    ev[col] = row.get(col, ".")
    # keep the first sample that saw each event (stable order)
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
                ev = events[eid]
                counts = [ev["counts"].get(s, 0) for s in args.sample_names]
                fh.write(eid + "\t" + "\t".join(str(c) for c in counts) + "\n")

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

    print(f"{len(order)} unique events across {len(args.sample_names)} samples "
          f"-> {args.out_events}")


if __name__ == "__main__":
    main()
