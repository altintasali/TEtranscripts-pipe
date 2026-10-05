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
                   STAR's per-sample read evidence is aggregated, not taken
                   from one sample: total_reads sums unique_reads,
                   multi_reads sums multi_reads, and overhang is the
                   maximum across samples -- the mapping-quality inputs for
                   candidates.tsv.gz's sj_unique_fraction / sj_max_overhang.
  counts_matrix.tsv  event_id x sample matrix of STAR's own unique_reads
                   count for that junction (0 where a sample never had this
                   junction). Written when --out-counts is given. Rows are
                   restricted to gene<->TE events (direction gene_to_te /
                   te_to_gene) -- see "Why the matrix is TE-restricted" below.
  cpm_matrix.tsv   the same matrix, each column divided by that sample's
                   column total x 1e6 (CPM, not TPM -- a splice junction has
                   no meaningful "length" to normalize by, same rationale as
                   chimera_chimeric_reads_counts.py). Written when --out-cpm is given.
                   Per-sample totals for the CPM denominator are summed over
                   ALL events (every direction), not just the TE-restricted
                   rows this file reports -- CPM is meant to read as "per
                   million reads in this sample's whole splice-junction
                   library", the standard library-size normalization; if the
                   denominator were also TE-restricted, CPM would instead
                   measure "share of this sample's TE-junction reads", which
                   answers a different question and would silently break
                   comparability with chimera_chimeric_reads_counts.py's CPM.
  te-gene-junctions.tsv  the all_events catalog filtered to gene<->TE events
                   (direction gene_to_te / te_to_gene), written when
                   --out-te-events is given. Same row set as counts_matrix/
                   cpm_matrix now cover.

Why the matrix is TE-restricted (measured 2026): counts_matrix/cpm_matrix
used to include every direction class STAR's SJ.out.tab produces --
gene_to_gene and other account for the large majority of junctions in a real
cohort (across 12 real mouse RNA-seq samples, gene<->TE events were ~4% of
all classified junctions: 76,430 of ~2M). Measured on real data (12 samples,
full genome annotation) and on a synthetic cohort built to match this
module's own worst-case estimate (~1M unique junctions, 80 samples): the
merge step itself stayed within its resources.yaml budget in both cases
(peak RSS 963 MB / 3.5 GB measured vs. 8 GB budgeted), but the downstream
DESeq2 vst() QC transform (sample_qc.R --transform) run on the unrestricted
matrix failed outright ("every gene contains at least one zero, cannot
compute log geometric means") and silently fell back to log2(counts + 1) --
a strictly worse normalization for the PCA/sample-distance QC view than a
real vst. Restricting the matrix to the ~13x-smaller gene<->TE subset (real
data: 23,910 vs 323,919 unique events at the same 12 samples) let vst()
actually succeed (6,399 features after the existing min_samples_present/
min_total_counts floor, 710 MB, 11s) instead of forcibly degrading to log2.
So the fix here is motivated by QC quality and by not holding ~13x more
matrix than the screen's own biology needs, not by an actual OOM risk under
the current resource budget -- measurement did not reproduce the OOM this
module was originally worried about.

Nothing is filtered here beyond what classify_chimera_splice_junctions.py itself already
applied (--min-unique-reads / --require-canonical) -- annotate-first, same as
the other two screens. The all_events catalog (--out-events) is NOT
TE-restricted: it stays the full annotate-first catalog across every
direction class, same as before, since holding it was measured cheap (see
above) and it is the natural place to look up a gene_to_gene/other junction
that MultiQC's Chimera section links to.
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
    "te_class", "chimera_type", "te_initiated_detail", "antisense_flag",
    "library_strand", "transcript_strand", "gene_strand_match",
    "gtf_annotated_intron",
]


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


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
            ev = events.setdefault(eid, {"sample": sample, "counts": {},
                                         "multi": 0, "max_overhang": 0})
            ev["counts"][sample] = int(row["unique_reads"])
            ev["multi"] += _int(row.get("multi_reads"))
            ev["max_overhang"] = max(ev["max_overhang"], _int(row.get("overhang")))
            if ev["sample"] == sample:
                for col in ANNOTATION_COLUMNS:
                    ev[col] = row.get(col, ".")
    for ev in events.values():
        ev["multi_reads"] = ev["multi"]
        ev["overhang"] = ev["max_overhang"]
    order = sorted(events, key=lambda e: (events[e]["sample"], e))
    te_order = [
        eid for eid in order
        if events[eid].get("direction") in ("gene_to_te", "te_to_gene")
    ]

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
            for eid in te_order:
                counts = events[eid]["counts"]
                vals = [counts.get(s, 0) for s in args.sample_names]
                fh.write(eid + "\t" + "\t".join(str(c) for c in vals) + "\n")

    if args.out_cpm:
        # Denominator uses ALL events (every direction), not just te_order --
        # see the module docstring's "Why the matrix is TE-restricted" note.
        totals = dict.fromkeys(args.sample_names, 0)
        for eid in order:
            counts = events[eid]["counts"]
            for s in args.sample_names:
                totals[s] += counts.get(s, 0)
        with open_write(args.out_cpm) as fh:
            fh.write("event_id\t" + "\t".join(args.sample_names) + "\n")
            for eid in te_order:
                counts = events[eid]["counts"]
                cpm = [
                    0 if totals[s] == 0 else counts.get(s, 0) / totals[s] * 1e6
                    for s in args.sample_names
                ]
                fh.write(eid + "\t" + "\t".join(f"{v:.3f}" for v in cpm) + "\n")

    if args.out_te_events:
        with open_write(args.out_te_events) as fh:
            fh.write("\t".join(ANNOTATION_COLUMNS + ["n_samples", "total_reads"]) + "\n")
            for eid in te_order:
                ev = events[eid]
                row = [ev.get(c, ".") for c in ANNOTATION_COLUMNS]
                row += [len(ev["counts"]), sum(ev["counts"].values())]
                fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{len(order)} unique junctions across {len(args.sample_names)} samples "
          f"-> {args.out_events}")


if __name__ == "__main__":
    main()
