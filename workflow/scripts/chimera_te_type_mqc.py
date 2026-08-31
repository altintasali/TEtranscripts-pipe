#!/usr/bin/env python3
"""TE type for the READ-evidence screen, per sample.

The assembly screen has its own TE-type view already
(chimera_assembly_summary_mqc.py's class chart); this is the matching one for
the read screen, so both screens answer the same question in the same shape.

The two are shown SEPARATELY on purpose. Both emit te_initiated /
te_exonized / te_terminated, so putting them in one plot invites reading
agreement as corroboration -- and they are not the same measurement:

  Read evidence (classify_chimera_junctions.py) classifies by GENOMIC POSITION:
  where the TE sits relative to the gene body, oriented by gene strand.
    te_initiated   TE lies upstream of the gene
    te_terminated  TE lies downstream
    te_exonized    TE overlaps / sits inside the gene span
  A trans event (two chromosomes) gets no class at all -- comparing
  coordinates across chromosomes is meaningless.

  Transcript evidence (classify_chimera_assembly.py) classifies by TRANSCRIPT
  STRUCTURE: which exon of the assembled transcript overlaps a TE (first /
  last / internal), and has two classes the read screen cannot produce.

So a TE inside a gene's intron that becomes the transcript's first exon is
te_initiated to the assembly screen and te_exonized to the read screen.
Neither is wrong. Each section says so.

Per sample, because the read screen has a per-sample number and this is where
an outlier shows up -- the removed per-sample status grid was the only other
place that view existed.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

CLASSES = ["te_initiated", "te_exonized", "te_terminated"]

CROSS_SCREEN_NOTE = (
    "<br><br><strong>The assembly screen uses these same words for a "
    "different measurement.</strong> Here a class is decided by <em>genomic "
    "position</em> &mdash; where the TE sits relative to the gene body: "
    "<code>te_initiated</code> upstream of the gene, "
    "<code>te_terminated</code> downstream, <code>te_exonized</code> inside "
    "it. In <strong>Assembly - TE type</strong> it is decided by "
    "<em>transcript structure</em> &mdash; whether the TE hits the first, "
    "last or an internal exon of the assembled transcript. A TE in a gene's "
    "intron that becomes the transcript's first exon is "
    "<code>te_initiated</code> there and <code>te_exonized</code> here, and "
    "neither is wrong. Do not read agreement between the two as "
    "corroboration."
)


def load_metrics(path):
    out = {}
    with open_read(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    return out


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-tables", required=True, nargs="+")
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.qc_tables) != len(args.samples):
        sys.exit("error: --qc-tables and --samples must have equal length")

    data = {}
    for path, sample in zip(args.qc_tables, args.samples):
        m = load_metrics(path)
        data[sample] = {c: _int(m.get(f"chimera_type_{c}")) for c in CLASSES}

    doc = {
        "id": "chimera_reads_te_type",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Reads - TE type",
        "description": (
            "Per-sample chimeric junctions by TE type, from STAR's chimeric "
            "reads. Only gene-TE junctions on one chromosome get a type: a "
            "trans event has no meaningful position relative to the gene."
            + CROSS_SCREEN_NOTE
        ),
    }
    if any(v for row in data.values() for v in row.values()):
        doc.update({
            "plot_type": "bargraph",
            "pconfig": {
                "id": "chimera_reads_te_type_plot",
                "title": "Reads: chimeric junctions by TE type",
                "ylab": "junctions",
                # cpswitch gives the share-of-sample view, which is the right
                # percentage for a composition. One control, not two.
                "cpswitch": True,
                "cpswitch_counts_label": "Junction counts",
                "cpswitch_percent_label": "% of the sample's typed junctions",
                "tt_decimals": 0,
            },
            "categories": CLASSES,
            "data": data,
        })
    else:
        # an all-zero bargraph makes MultiQC exit non-zero and write no report
        doc.update({
            "plot_type": "html",
            "data": ("<p>No chimeric junction in this run received a TE type. "
                     "This is not an error: it means no gene-TE junction had "
                     "both partners on one chromosome.</p>"),
        })

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    totals = {c: sum(row[c] for row in data.values()) for c in CLASSES}
    print(f"reads TE type: {len(data)} samples, {totals} -> {args.out}")


if __name__ == "__main__":
    main()
