#!/usr/bin/env python3
"""Merge the per-sample junction QC tables (junction_qc.py) into MultiQC
custom-content bar plots:
  junction_qc_mqc.json      per-sample chimeric-junction composition by
                            direction (counts and % of total junctions).
  te_chimeras_mqc.json      the gene<->TE subset by direction (gene_to_te /
                            te_to_gene, counts and % of total junctions) --
                            the TE-chimeras view, written when --out-te-chimeras
                            is given.

Reads every results/chimera/qc/{sample}_junction_qc.tsv (metric/value pairs)
and writes the JSONs above. MultiQC renders them inside multiqc_report.html in
the custom_content module (ordered by multiqc_config.yaml), with the two
datasets switchable via each plot's cpswitch control.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _io import open_read, open_write

DIRECTIONS = [
    "gene_to_te", "te_to_gene", "gene_to_gene", "te_to_te",
    "gene_to_other", "other_to_gene", "te_to_other", "other",
]


def load_metrics(path):
    """metric/value pairs from a junction_qc.tsv."""
    metrics = {}
    with open_read(path) as fh:
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
    ap.add_argument(
        "--out-te-chimeras", required=False,
        help="Optional second output: the gene<->TE subset (direction "
        "gene_to_te / te_to_gene) as its own bar plot.",
    )
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
            "data_labels": [
                {"name": "Junction counts"},
                {"name": "% of total junctions"},
            ],
        },
        "data": [counts, pct],
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"junction QC barplot ({len(samples)} samples) -> {args.out}")

    if args.out_te_chimeras:
        te_dirs = ["gene_to_te", "te_to_gene"]
        te_counts = {}
        metrics_by_sample = {}
        for path, sample in zip(args.tables, args.samples):
            metrics = load_metrics(path)
            metrics_by_sample[sample] = metrics
            per_sample = {}
            for d in te_dirs:
                try:
                    per_sample[d] = int(float(metrics.get(f"direction_{d}", 0)))
                except (TypeError, ValueError):
                    per_sample[d] = 0
            te_counts[sample] = per_sample

        te_pct = {}
        for sample in samples:
            try:
                total = int(float(metrics_by_sample[sample].get("events_total", 0)))
            except (TypeError, ValueError):
                total = 0
            if total <= 0:
                te_pct[sample] = {d: 0.0 for d in te_dirs}
            else:
                te_pct[sample] = {
                    d: round(te_counts[sample][d] * 100.0 / total, 1)
                    for d in te_dirs
                }

        te_doc = {
            "id": "chimera_te_chimeras",
            "section_name": "TE chimeras",
            "description": (
                "Per-sample gene\u2194TE chimeric junctions (direction "
                "gene_to_te / te_to_gene), as counts and % of total junctions."
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_te_chimeras_plot",
                "title": "Gene-TE chimeras by direction",
                "ylab": "junctions",
                "cpswitch": True,
                "data_labels": [
                    {"name": "TE-chimera junctions"},
                    {"name": "% of total junctions"},
                ],
            },
            "data": [te_counts, te_pct],
        }

        os.makedirs(os.path.dirname(args.out_te_chimeras), exist_ok=True)
        with open_write(args.out_te_chimeras) as fh:
            json.dump(te_doc, fh, indent=2)
            fh.write("\n")
        print(f"TE-chimeras barplot ({len(samples)} samples) -> {args.out_te_chimeras}")


if __name__ == "__main__":
    main()
