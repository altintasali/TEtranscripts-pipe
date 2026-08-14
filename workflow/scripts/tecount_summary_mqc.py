#!/usr/bin/env python3
"""Per-sample TEcounts summary stats for the MultiQC report (custom content):
gene-vs-TE assignment and TE class composition.

TEcount writes one {sample}.cntTable per sample: a two-column table
(`gene/TE` feature key + read count) holding gene ids and TE subfamily names
mixed in one column. TEcount builds TE keys as `gene_id:family_id:class_id`
from the TE GTF, so each row can be classified from the key alone: a key with
no colons is a gene; a `gene:family:class` key is a TE subfamily whose family
and class come from the key.

Emits two MultiQC bar-plot custom-content documents:
  tecount_assignment_mqc.json  genes vs TE subfamilies (counts + % of total)
  tecount_te_class_mqc.json    TE subfamily reads by class, LINE/SINE/LTR/
                               DNA/RC/unknown (counts + % of TE reads)

Unlike the sample-QC matrix (tecount_counts.py), this uses the RAW cntTables
(all features), so it is independent of tetranscripts.qc.feature_class.
"""
import argparse
import json
import os
import sys

TE_CLASSES = ["LINE", "SINE", "LTR", "DNA", "RC"]


def load_counts(path):
    """feature -> count from a raw TEcount cntTable (header + rows)."""
    counts = {}
    with open(path) as fh:
        fh.readline()  # header: gene/TE \t <bam path>
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                counts[parts[0]] = int(float(parts[1]))
            except ValueError:
                counts[parts[0]] = 0
    return counts


def classify(key):
    """Return ("gene", None, None) or ("te", family, class)."""
    parts = key.split(":")
    if len(parts) >= 2 and all(parts):
        family = parts[-2]
        cls = parts[-1] if len(parts) >= 3 else "unknown"
        if cls not in TE_CLASSES:
            cls = "unknown"
        return ("te", family, cls)
    return ("gene", None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--out-assignment", required=True)
    ap.add_argument("--out-class", required=True)
    args = ap.parse_args()

    if len(args.tables) != len(args.samples):
        sys.exit("error: --tables and --samples must have equal length")

    assign = {}   # sample -> {"genes": n, "TE subfamilies": n}
    classes = {}  # sample -> {class: n}
    for path, sample in zip(args.tables, args.samples):
        gene, te = 0, 0
        cls_counts = {}
        for key, n in load_counts(path).items():
            kind, _family, cls = classify(key)
            if kind == "gene":
                gene += n
            else:
                te += n
                cls_counts[cls] = cls_counts.get(cls, 0) + n
        assign[sample] = {"genes": gene, "TE subfamilies": te}
        classes[sample] = {c: cls_counts.get(c, 0) for c in TE_CLASSES}
        classes[sample]["unknown"] = cls_counts.get("unknown", 0)

    samples = sorted(assign)

    assign_pct = {}
    for s in samples:
        total = assign[s]["genes"] + assign[s]["TE subfamilies"]
        if total <= 0:
            assign_pct[s] = {k: 0.0 for k in assign[s]}
        else:
            assign_pct[s] = {
                k: round(v * 100.0 / total, 1)
                for k, v in assign[s].items()
            }

    class_pct = {}
    for s in samples:
        total = sum(classes[s].values())
        if total <= 0:
            class_pct[s] = {c: 0.0 for c in classes[s]}
        else:
            class_pct[s] = {
                c: round(v * 100.0 / total, 1)
                for c, v in classes[s].items()
            }

    assignment_doc = {
        "id": "tecount_assignment",
        "section_name": "TEcounts: gene vs TE assignment",
        "description": (
            "Per-sample read counts assigned by TEcount to genes vs TE "
            "subfamilies (counts and % of total assigned reads)."
        ),
        "plot_type": "bar",
        "pconfig": {
            "id": "tecount_assignment_plot",
            "title": "TEcounts assignment (genes vs TEs)",
            "ylab": "assigned reads",
            "cpswitch": True,
            "data_labels": [
                {"name": "Read counts"},
                {"name": "% of total assigned reads"},
            ],
        },
        "data": [assign, assign_pct],
    }
    class_doc = {
        "id": "tecount_te_class",
        "section_name": "TEcounts: TE class composition",
        "description": (
            "Per-sample TE-subfamily reads by repeat class (counts and % of "
            "TE reads); classes not in LINE/SINE/LTR/DNA/RC are grouped as "
            "unknown."
        ),
        "plot_type": "bar",
        "pconfig": {
            "id": "tecount_te_class_plot",
            "title": "TE class composition",
            "ylab": "TE reads",
            "cpswitch": True,
            "data_labels": [
                {"name": "TE read counts"},
                {"name": "% of TE reads"},
            ],
        },
        "data": [classes, class_pct],
    }

    for out, doc in ((args.out_assignment, assignment_doc),
                     (args.out_class, class_doc)):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")
        print(f"TEcounts summary barplot ({len(samples)} samples) -> {out}")


if __name__ == "__main__":
    main()
