#!/usr/bin/env python3
"""Per-sample TElocal summary stats for the MultiQC report (custom content):
gene-vs-TE assignment and TE class composition.

TElocal writes one {sample}.cntTable per sample: a two-column table
(`gene/TE` feature key + read count) holding gene ids and TE locus keys
mixed in one column. TElocal builds TE keys as
`transcript_id:gene_id:family_id:class_id` from the .locInd index, where
transcript_id embeds genomic coordinates (e.g.
`chr1:564318:564741(L1PA2:+):L1PA2:L1:LINE/L1`). Each row can be classified
from the key alone: a key with no colons is a gene; a key with >=2 colons
has family = second-to-last field, class = last field.

Emits two MultiQC bar-plot custom-content documents:
  telocal_assignment_mqc.json  genes vs TE loci (counts + % of total)
  telocal_te_class_mqc.json    TE locus reads by class, LINE/SINE/LTR/
                               DNA/RC/unknown (counts + % of TE reads)

Unlike the sample-QC matrix, this uses the RAW cntTables (all features).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

TE_CLASSES = ["LINE", "SINE", "LTR", "DNA", "RC"]


def load_counts(path):
    """feature -> count from a raw TElocal cntTable (header + rows)."""
    counts = {}
    with open_read(path) as fh:
        fh.readline()  # header: gene/TE \t <bam path>
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                counts[parts[0]] = int(parts[1])
            except ValueError:
                counts[parts[0]] = 0
    return counts


def classify(key):
    """Return ("gene", None, None) or ("te", family, class).

    TElocal TE keys have the form
    transcript_id:gene_id:family_id:class_id  (4+ colon-separated fields).
    The last two fields are always family and class, matching TEcount's
    convention. Gene keys have no colons.
    """
    parts = key.split(":")
    if len(parts) >= 3 and all(parts):
        family = parts[-2]
        cls = parts[-1] if len(parts) >= 2 else "unknown"
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

    assign = {}   # sample -> {"genes": n, "TE loci": n}
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
        assign[sample] = {"genes": gene, "TE loci": te}
        classes[sample] = {c: cls_counts.get(c, 0) for c in TE_CLASSES}
        classes[sample]["unknown"] = cls_counts.get("unknown", 0)

    samples = sorted(assign)

    assign_pct = {}
    for s in samples:
        total = assign[s]["genes"] + assign[s]["TE loci"]
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
        "id": "telocal_assignment",
        "parent_id": "telocal",
        "parent_name": "TElocal",
        "section_name": "Gene vs. TE Assignment",
        "description": (
            "Per-sample read counts assigned by TElocal to genes vs TE "
            "loci (counts and % of total assigned reads). "
            "<br><br><em>How to read this:</em> TElocal resolves each TE "
            "insertion separately, where TEcount pools all copies of a "
            "subfamily into a single row. A subfamily that looks uniformly "
            "expressed in TEcount is often one or two highly-expressed "
            "copies here, with the rest silent \u2014 which is the reason "
            "to run both."
        ),
        "plot_type": "bar",
        "pconfig": {
            "id": "telocal_assignment_plot",
            "title": "TElocal assignment (genes vs TEs)",
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
        "id": "telocal_te_class",
        "parent_id": "telocal",
        "parent_name": "TElocal",
        "section_name": "TE Class Composition",
        "description": (
            "Per-sample TE-locus reads by repeat class (counts and % of "
            "TE reads); classes not in LINE/SINE/LTR/DNA/RC are grouped as "
            "unknown."
        ),
        "plot_type": "bar",
        "pconfig": {
            "id": "telocal_te_class_plot",
            "title": "TElocal TE class composition",
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
        with open_write(out) as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")
        print(f"TElocal summary barplot ({len(samples)} samples) -> {out}")


if __name__ == "__main__":
    main()
