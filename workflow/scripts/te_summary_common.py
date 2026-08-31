#!/usr/bin/env python3
"""Shared implementation behind the TEcount and TElocal summary sections.

These two sections answer the same question at two resolutions -- how many
reads went to genes vs TEs, and how the TE reads split by repeat class -- from
the same two-column cntTable format. They were two near-identical scripts, 88
differing lines out of 358, and they had drifted in a way nobody would notice:

  TEcount:  counts[key] = int(float(value))
  TElocal:  counts[key] = int(value)          <- ValueError, caught, -> 0

TEcount's EM assigns fractional reads before rounding, so a cntTable can carry
"12.0". The TElocal path turned every such row into a silent ZERO rather than
failing, which is data loss inside a report that looks fine. One tolerant
loader lives here now.

What genuinely differs between the two stays in FLAVOURS: the key format (how
many colon-separated fields make a TE), the section ids and names, and the
prose. Everything else is shared.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

TE_CLASSES = ["LINE", "SINE", "LTR", "DNA", "RC"]


class Flavour:
    """Everything that differs between the TEcount and TElocal sections."""

    def __init__(self, tool, parent_id, unit_label, min_key_fields,
                 assignment_desc, class_desc, assignment_title, class_title):
        self.tool = tool
        # parent_id must match sample_qc.R's view id for this tool -- a
        # mismatch silently splits the group into two report sections.
        self.parent_id = parent_id
        self.unit_label = unit_label          # bar series name, e.g. "TE loci"
        self.min_key_fields = min_key_fields  # colon-separated fields => a TE
        self.assignment_desc = assignment_desc
        self.class_desc = class_desc
        self.assignment_title = assignment_title
        self.class_title = class_title


def load_counts(path):
    """feature -> count from a raw cntTable (header + rows).

    int(float(...)) deliberately: counts can arrive as "12.0" from the EM
    step. A row that is genuinely unparseable still falls back to 0, but a
    float-formatted integer must not.
    """
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
                counts[parts[0]] = int(float(parts[1]))
            except ValueError:
                counts[parts[0]] = 0
    return counts


def classify(key, min_fields):
    """Return ("gene", None, None) or ("te", family, class).

    Both tools put family and class in the LAST TWO colon-separated fields;
    they differ only in how many fields a TE key has. TEcount writes
    gene_id:family_id:class_id, TElocal transcript_id:gene_id:family_id:
    class_id (where transcript_id itself embeds coordinates and colons).
    A key with too few fields is a gene id.
    """
    parts = key.split(":")
    if len(parts) >= min_fields and all(parts):
        family = parts[-2]
        # A 2-field key has no class field of its own -- the last field IS
        # the family -- so the class is unknown rather than a copy of it.
        # (Only reachable for TEcount, whose min_fields is 2.)
        cls = parts[-1] if len(parts) >= 3 else "unknown"
        if cls not in TE_CLASSES:
            cls = "unknown"
        return ("te", family, cls)
    return ("gene", None, None)


def summarise(flavour, tables, samples):
    """(assignment counts, class counts) per sample."""
    assign, classes = {}, {}
    for path, sample in zip(tables, samples):
        gene, te = 0, 0
        cls_counts = {}
        for key, n in load_counts(path).items():
            kind, _family, cls = classify(key, flavour.min_key_fields)
            if kind == "gene":
                gene += n
            else:
                te += n
                cls_counts[cls] = cls_counts.get(cls, 0) + n
        assign[sample] = {"genes": gene, flavour.unit_label: te}
        classes[sample] = {c: cls_counts.get(c, 0) for c in TE_CLASSES}
        classes[sample]["unknown"] = cls_counts.get("unknown", 0)
    return assign, classes


def _as_percent(per_sample):
    """Row-normalise to percentages; an all-zero sample stays all-zero."""
    out = {}
    for sample, values in per_sample.items():
        total = sum(values.values())
        if total <= 0:
            out[sample] = dict.fromkeys(values, 0.0)
        else:
            out[sample] = {k: round(v * 100.0 / total, 1)
                           for k, v in values.items()}
    return out


def build_docs(flavour, assign, classes):
    """The two MultiQC custom-content documents for this flavour."""
    assignment_doc = {
        "id": f"{flavour.parent_id}_assignment",
        "parent_id": flavour.parent_id,
        "parent_name": flavour.tool,
        "section_name": "Gene vs. TE Assignment",
        "description": flavour.assignment_desc,
        "plot_type": "bar",
        "pconfig": {
            "id": f"{flavour.parent_id}_assignment_plot",
            "title": flavour.assignment_title,
            "ylab": "assigned reads",
            "cpswitch": True,
            "data_labels": [
                {"name": "Read counts", "format": "{:,.0f}"},
                {"name": "% of total assigned reads", "format": "{:,.1f}"},
            ],
        },
        "data": [assign, _as_percent(assign)],
    }
    class_doc = {
        "id": f"{flavour.parent_id}_te_class",
        "parent_id": flavour.parent_id,
        "parent_name": flavour.tool,
        "section_name": "TE Class Composition",
        "description": flavour.class_desc,
        "plot_type": "bar",
        "pconfig": {
            "id": f"{flavour.parent_id}_te_class_plot",
            "title": flavour.class_title,
            "ylab": "TE reads",
            "cpswitch": True,
            "data_labels": [
                {"name": "TE read counts", "format": "{:,.0f}"},
                {"name": "% of TE reads", "format": "{:,.1f}"},
            ],
        },
        "data": [classes, _as_percent(classes)],
    }
    return assignment_doc, class_doc


def run(flavour, args):
    """Shared entry point: validate, summarise, write both documents."""
    if len(args.tables) != len(args.samples):
        sys.exit("error: --tables and --samples must have equal length")

    assign, classes = summarise(flavour, args.tables, args.samples)
    assignment_doc, class_doc = build_docs(flavour, assign, classes)

    for out, doc in ((args.out_assignment, assignment_doc),
                     (args.out_class, class_doc)):
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        with open_write(out) as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")
        print(f"{flavour.tool} summary barplot "
              f"({len(sorted(assign))} samples) -> {out}")
