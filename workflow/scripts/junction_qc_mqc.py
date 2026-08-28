#!/usr/bin/env python3
"""Merge the per-sample junction QC tables (junction_qc.py) into MultiQC
custom-content bar plots:
  junction_qc_mqc.json      per-sample chimeric-junction composition by
                            direction (counts and % of total junctions).
  te_gene_chimeras_mqc.json the gene<->TE subset by direction (gene_to_te /
                            te_to_gene, counts and % of total junctions) --
                            the gene-TE chimeras view, written when
                            --out-te-gene-chimeras is given.

Reads every results/chimera_junction/qc/{sample}_junction_qc.tsv (metric/value pairs)
and writes the JSONs above. MultiQC renders them inside multiqc_report.html in
the custom_content module (ordered by multiqc_config.yaml), with the two
datasets switchable via each plot's cpswitch control.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

DIRECTIONS = [
    "gene_to_te", "te_to_gene", "gene_to_gene", "te_to_te",
    "gene_to_other", "other_to_gene", "te_to_other", "other_to_te", "other",
]


def _finalise(doc, datasets, empty_html):
    """Attach bar datasets, or fall back to html when everything is zero.

    MultiQC's bargraph raises "No datasets to plot" on an all-zero document
    and takes the whole report down with it. A run with no junctions -- or,
    for the canonical plot, no junction carrying a splice motif -- is a
    perfectly normal outcome (small genome, shallow data, synthetic reads),
    so it must be documented rather than fatal.
    """
    if any(v for row in datasets[0].values() for v in row.values()):
        doc["plot_type"] = "bar"
        doc["data"] = datasets
    else:
        doc["plot_type"] = "html"
        doc["data"] = empty_html
        doc.pop("pconfig", None)
    return doc


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
        "--out-canonical", required=False,
        help="Optional output: canonical (splice-motif) rate per direction "
        "-- the main signal-vs-artifact discriminator.",
    )
    ap.add_argument(
        "--out-te-gene-chimeras", required=False,
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
        "parent_id": "chimera_reads",
        "parent_name": "Gene-TE chimeras: read evidence",
        "section_name": "Events by direction",
        "description": (
            "Per-sample composition of annotated chimeric junctions by "
            "direction (the gene\u2194TE event classes first), as counts and "
            "% of total junctions. "
            "<br><br><em>How to read this:</em> a junction is classified by "
            "what its two breakpoints overlap. Only <code>gene_to_te</code> "
            "and <code>te_to_gene</code> are gene\u2013TE chimeras. The other "
            "classes are not merely leftovers \u2014 STAR calls a junction "
            "chimeric on alignment geometry alone, without reading any "
            "annotation, so they also collect circRNA back-splices, "
            "read-through transcripts and PCR/ligation chimeras. Nothing is "
            "filtered here; see the canonical-rate plot below before "
            "treating any class as real."
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
    }
    _finalise(doc, [counts, pct],
              "<p>No chimeric junctions were annotated in this run, so there "
              "is nothing to plot. This is not an error \u2014 STAR found no "
              "chimeric reads, or none survived annotation.</p>")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"junction QC barplot ({len(samples)} samples) -> {args.out}")

    if args.out_canonical:
        canon_pct, canon_n = {}, {}
        for path, sample in zip(args.tables, args.samples):
            m = load_metrics(path)
            pct_row, n_row = {}, {}
            for d in DIRECTIONS:
                try:
                    tot = int(float(m.get(f"direction_{d}", 0)))
                    hit = int(float(m.get(f"canonical_{d}", 0)))
                except (TypeError, ValueError):
                    tot = hit = 0
                pct_row[d] = round(100.0 * hit / tot, 1) if tot else 0.0
                n_row[d] = hit
            canon_pct[sample], canon_n[sample] = pct_row, n_row

        canon_doc = {
            "id": "chimera_canonical_rate",
            "parent_id": "chimera_reads",
            "parent_name": "Gene-TE chimeras: read evidence",
            "section_name": "Splice-motif rate by direction",
            "description": (
                "Share of each direction's junctions for which STAR reported a "
                "recognised splice motif (GT/AG, GC/AG, AT/AC and reverse "
                "complements). "
                "<br><br><em>Why this matters:</em> genuine spliced introns are "
                "very nearly 100% canonical, so a junction with no motif is "
                "most likely template switching, a ligation/PCR chimera or a "
                "mismapping \u2014 not a transcript. Chimeric junctions are "
                "overwhelmingly motif-less in practice, so <strong>read this "
                "as enrichment, not as an absolute</strong>. "
                "<br><br><strong>Compare within a donor group, not against "
                "<code>other</code>.</strong> The donor side alone shifts the "
                "rate a lot \u2014 on real data every <code>gene_to_*</code> "
                "class sits near 8-10% while every <code>te_to_*</code> class "
                "sits at 3-7% \u2014 so measuring <code>gene_to_te</code> "
                "against <code>other</code> credits it for the donor being a "
                "gene at all. The honest comparisons are "
                "<code>gene_to_te</code> vs <code>gene_to_gene</code> / "
                "<code>gene_to_other</code>, and <code>te_to_gene</code> vs "
                "<code>te_to_te</code> / <code>te_to_other</code>. Both "
                "gene\u2013TE directions beating their own donor group is "
                "the evidence the classification captures something real; "
                "the margin is modest (well under 2x), not several-fold. "
                "The canonical subset of the gene\u2013TE classes is the "
                "sensible working set."
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_canonical_rate_plot",
                "title": "Canonical (splice-motif) rate by direction",
                "ylab": "% canonical",
                "cpswitch": False,
                "data_labels": [
                    {"name": "% canonical"},
                    {"name": "canonical junctions"},
                ],
            },
        }
        _finalise(canon_doc, [canon_pct, canon_n],
                  "<p>No chimeric junction in this run carried a recognised "
                  "splice motif, so there is no canonical rate to plot. That "
                  "is expected on synthetic or very shallow data; on a real "
                  "library it means the chimeric calls are dominated by "
                  "motif-less artifacts.</p>")
        os.makedirs(os.path.dirname(args.out_canonical), exist_ok=True)
        with open_write(args.out_canonical) as fh:
            json.dump(canon_doc, fh, indent=2)
            fh.write("\n")
        print(f"canonical-rate barplot ({len(samples)} samples) -> {args.out_canonical}")

    if args.out_te_gene_chimeras:
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
            "id": "chimera_te_gene_chimeras",
            "parent_id": "chimera_reads",
            "parent_name": "Gene-TE chimeras: read evidence",
            "section_name": "Gene-TE subset",
            "description": (
                "Per-sample gene\u2194TE chimeric junctions (direction "
                "gene_to_te / te_to_gene), as counts and % of total junctions. "
                "<br><br><em>How to read this:</em> these are candidate "
                "gene\u2013TE chimeras, annotated but <strong>not</strong> "
                "filtered \u2014 no read-count, replicate or splice-motif "
                "cutoff has been applied. Apply your own before treating a "
                "call as confident, and check the canonical-rate plot: in "
                "practice only a minority of chimeric junctions carry a "
                "splice motif at all."
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_te_gene_chimeras_plot",
                "title": "Gene-TE chimeras by direction",
                "ylab": "junctions",
                "cpswitch": True,
                "data_labels": [
                    {"name": "Gene-TE chimeric junctions"},
                    {"name": "% of total junctions"},
                ],
            },
        }
        _finalise(te_doc, [te_counts, te_pct],
                  "<p>No gene\u2013TE chimeric junctions were found in this "
                  "run, so there is nothing to plot.</p>")

        os.makedirs(os.path.dirname(args.out_te_gene_chimeras), exist_ok=True)
        with open_write(args.out_te_gene_chimeras) as fh:
            json.dump(te_doc, fh, indent=2)
            fh.write("\n")
        print(f"Gene-TE chimeras barplot ({len(samples)} samples) -> {args.out_te_gene_chimeras}")


if __name__ == "__main__":
    main()
