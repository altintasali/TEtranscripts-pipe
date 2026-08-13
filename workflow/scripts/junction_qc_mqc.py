#!/usr/bin/env python3
"""Merge the per-sample junction QC tables (junction_qc.py) into one MultiQC
custom-content bar plot: per-sample chimeric-junction composition by direction
(counts and % of total junctions).

Reads every results/chimera/qc/{sample}_junction_qc.tsv (metric/value pairs)
and writes results/chimera/qc/junction_qc_mqc.json. MultiQC renders it inside
multiqc_report.html in the custom_content module (ordered by
multiqc_config.yaml), with the two datasets switchable via the plot's
cpswitch control.
"""
import argparse
import json
import os
import sys

DIRECTIONS = [
    "gene_to_te", "te_to_gene", "gene_to_gene", "te_to_te",
    "gene_to_other", "other_to_gene", "te_to_other", "other",
]


def load_metrics(path):
    """metric/value pairs from a junction_qc.tsv."""
    metrics = {}
    with open(path) as fh:
        fh.readline()  # header: metric \t value
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                metrics[parts[0]] = parts[1]
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.tables) != len(args.samples):
        sys.exit("error: --tables and --samples must have equal length")

    counts = {}
    for path, sample in zip(args.tables, args.samples):
        metrics = load_metrics(path)
        per_sample = {}
        for d in DIRECTIONS:
            try:
                per_sample[d] = int(float(metrics.get(f"direction_{d}", 0)))
            except (TypeError, ValueError):
                per_sample[d] = 0
        counts[sample] = per_sample

    samples = sorted(counts)
    pct = {}
    for sample in samples:
        total = sum(counts[sample].values())
        if total <= 0:
            pct[sample] = {d: 0.0 for d in DIRECTIONS}
        else:
            pct[sample] = {
                d: round(counts[sample][d] * 100.0 / total, 1)
                for d in DIRECTIONS
            }

    doc = {
        "id": "chimera_junction_qc",
        "section_name": "Chimera junction QC",
        "description": (
            "Per-sample composition of annotated chimeric junctions by "
            "direction (the gene\u2194TE event classes first), as counts and "
            "% of total junctions."
        ),
        "plot_type": "bar",
        "pconfig": {
            "id": "chimera_junction_qc_plot",
            "title": "Chimera junctions by direction",
            "ylab": "junctions",
            "cpswitch": True,
        },
        "samples": samples,
        "datasets": [
            {"name": "Junction counts", "data": counts},
            {"name": "% of total junctions", "data": pct},
        ],
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"junction QC barplot ({len(samples)} samples) -> {args.out}")


if __name__ == "__main__":
    main()
