#!/usr/bin/env python3
"""Per-sample evidence status grid, in the shape of FastQC's Status Checks:
samples down the side, evidence layers across the top, one traffic-light cell
each.

The other two evidence heatmaps summarise the cohort -- they answer "which
evidence types agree" and "which candidates lead". Neither answers the
question you ask first when a run finishes: **is one of my samples
different?** A sample contributing a tenth of the gene-TE junctions its
replicates do, or the only one with a strandedness mismatch, is invisible in
any cohort-level plot because the cohort absorbs it.

Thresholds are RELATIVE TO THE COHORT, not absolute. Absolute cut-offs are
meaningless here: chimeric-junction yield varies by orders of magnitude with
depth, library prep and genome, so "500 gene-TE junctions" is unremarkable in
one experiment and alarming in another. What is always meaningful is one
sample departing from its peers, so each layer is scored against the cohort
median:

  pass   within 2x of the median (or the median is 0 -- nothing to compare)
  warn   2-5x away
  fail   more than 5x away, or zero where the cohort median is not

With fewer than 3 samples there is no usable median, so every cell is
reported as pass and the description says so rather than inventing a verdict
from n=2.

Strandedness is the exception and is taken as-is from the strandedness check:
a MISMATCH is a fail regardless of what the other samples did, because it is
a statement about correctness rather than about yield.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read

PASS, WARN, FAIL = 1.0, 0.5, 0.0

# (label, how to derive the per-sample number from the junction QC metrics)
LAYERS = [
    ("Chimeric junctions", lambda m: m("events_total")),
    ("Gene-TE junctions", lambda m: m("events_gene_te")),
    ("Splice motif", lambda m: m("canonical_gene_to_te") + m("canonical_te_to_gene")),
    ("Unambiguous direction",
     lambda m: m("events_gene_te") - m("ambiguous_direction")),
    ("TE-initiated", lambda m: m("chimera_type_te_initiated")),
    ("TE-exonized", lambda m: m("chimera_type_te_exonized")),
    ("TE-terminated", lambda m: m("chimera_type_te_terminated")),
    ("Gene strand match", lambda m: m("strand_match_yes")),
]


def load_metrics(path):
    out = {}
    with open_read(path) as fh:
        fh.readline()  # metric \t value
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    return out


def getter(metrics):
    def g(key):
        try:
            return float(metrics.get(key, 0) or 0)
        except ValueError:
            return 0.0
    return g


def median(values):
    v = sorted(values)
    n = len(v)
    if not n:
        return 0.0
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def status(value, med):
    """Traffic light for one cell, scored against the cohort median."""
    if med <= 0:
        # Nothing to compare against -- every sample is equally empty. Saying
        # "fail" here would flag the whole column red for a layer that simply
        # does not apply to this experiment.
        return PASS
    if value <= 0:
        return FAIL
    ratio = value / med if value >= med else med / value
    if ratio <= 2:
        return PASS
    return WARN if ratio <= 5 else FAIL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-tables", required=True, nargs="+")
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--strandedness", default=None,
                    help="strandedness_check_mqc.json, for the extra column")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    samples = list(args.samples)
    metrics = {s: getter(load_metrics(p))
               for s, p in zip(samples, args.qc_tables)}

    labels = [lbl for lbl, _ in LAYERS]
    raw = [[fn(metrics[s]) for _, fn in LAYERS] for s in samples]
    meds = [median([raw[r][c] for r in range(len(samples))])
            for c in range(len(LAYERS))]
    enough = len(samples) >= 3
    grid = [[status(raw[r][c], meds[c]) if enough else PASS
             for c in range(len(LAYERS))] for r in range(len(samples))]

    strand_note = ""
    if args.strandedness and os.path.exists(args.strandedness):
        with open(args.strandedness) as fh:
            doc = json.load(fh)
        rows = doc.get("data", {})
        labels.append("Strandedness")
        n_mismatch = 0
        for r, s in enumerate(samples):
            st = rows.get(s, {}).get("status", "")
            if st == "MISMATCH":
                grid[r].append(FAIL)
                n_mismatch += 1
            elif st == "OK":
                grid[r].append(PASS)
            else:  # auto-detected: nothing was declared to disagree with
                grid[r].append(WARN)
        strand_note = (
            f" <strong>{n_mismatch} sample(s) declare a strandedness the data "
            "disagrees with</strong> -- see the Strandedness check section."
            if n_mismatch else
            " Strandedness agrees with the sample sheet everywhere."
        )

    scale = (
        "Thresholds are relative to the cohort median for each column, since "
        "absolute chimeric-junction yield varies by orders of magnitude with "
        "depth and library prep: <strong>green</strong> within 2x of the "
        "median, <strong>amber</strong> 2-5x away, <strong>red</strong> "
        "beyond 5x or zero where the cohort is not."
        if enough else
        "<strong>Fewer than 3 samples, so there is no usable cohort median "
        "and every cell is reported as pass.</strong> This grid detects "
        "outliers among replicates; it cannot do that for n &lt; 3."
    )
    doc = {
        "id": "sample_evidence_status",
        "parent_id": "evidence_status",
        "parent_name": "Evidence status",
        "section_name": "Per-sample status",
        "description": (
            "Which evidence layers each sample contributed to, in the shape "
            "of FastQC's status checks. Read it across a row to see what one "
            "sample produced, and down a column to spot the sample that "
            f"disagrees with its replicates. {scale}{strand_note}"
        ),
        "plot_type": "heatmap",
        "pconfig": {
            "id": "sample_evidence_status_plot",
            "title": "Evidence layers per sample",
            "min": 0, "max": 1,
            "colstops": [[0, "#d9534f"], [0.5, "#f0ad4e"], [1, "#5cb85c"]],
            "reverseColors": False,
        },
        "xcats": labels,
        "ycats": samples,
        "data": grid,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    n_flag = sum(1 for r in grid for v in r if v < PASS)
    print(f"sample evidence status: {len(samples)} samples x {len(labels)} "
          f"layers, {n_flag} non-pass cell(s) -> {args.out}")


if __name__ == "__main__":
    main()
