#!/usr/bin/env python3
"""TE type for the READ-evidence screen, per sample.

The assembly screen has its own TE-type view already
(chimera_assembly_summary_mqc.py's class chart); this is the matching one for
the read screen, so both screens answer the same question in the same shape.

The two are shown SEPARATELY on purpose. Both emit te_initiated /
te_exonized / te_terminated, so putting them in one plot invites reading
agreement as corroboration -- and they are not the same measurement:

  Read evidence (classify_chimera_chimeric_reads.py) classifies by the
  junction's DIRECTION plus the gene's own annotated exon structure --
  matching the splice-junctions screen's own typing exactly (see that
  script's module docstring for the full reasoning):
    te_initiated   TE-to-gene junction (TE upstream of the exon it splices
                   into, whatever its genomic coordinate relative to the
                   gene's overall span)
    te_terminated  gene-to-TE junction with no other annotated exon further
                   downstream than the donor exon
    te_exonized    gene-to-TE junction with an annotated exon further
                   downstream (the TE sits inside a region the gene's
                   transcript is annotated to continue past)
  A trans event (two chromosomes) gets no class at all -- comparing
  coordinates across chromosomes is meaningless.

  Transcript evidence (classify_chimera_assembly.py) classifies by TRANSCRIPT
  STRUCTURE: which exon of the assembled transcript overlaps a TE (first /
  last / internal), and has two classes the read screen cannot produce.

So a TE inside a gene's intron that becomes the transcript's first exon
without ever being observed splicing directly into that exon on a single
chimeric read is te_initiated to the assembly screen (from the assembled
transcript's own TSS) but only te_exonized here if this screen never saw a
te_to_gene junction landing on it. Each section says so.

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

# antisense_to_gene: a gene<->TE event on the strand opposite its gene
# (stranded libraries only) -- see chimera_exon_context.ANTISENSE_TO_GENE.
CLASSES = ["te_initiated", "te_exonized", "te_terminated", "antisense_to_gene"]

CROSS_SCREEN_NOTE = (
    "<br><br><strong>The assembly screen uses these same words for a "
    "different measurement.</strong> Here a class is decided by the "
    "junction's <em>direction</em> plus the gene's own annotated exon "
    "structure (matching the splice-junctions screen exactly): "
    "<code>te_initiated</code> for a TE-to-gene junction, "
    "<code>te_terminated</code> / <code>te_exonized</code> for a "
    "gene-to-TE junction depending on whether an annotated exon lies "
    "further downstream of the donor exon. In <strong>Assembly - TE "
    "type</strong> it is decided by <em>transcript structure</em> "
    "&mdash; whether the TE hits the first, last or an internal exon of "
    "the assembled transcript. A TE in a gene's intron that becomes the "
    "transcript's first exon is <code>te_initiated</code> there, and "
    "<code>te_initiated</code> here only if a chimeric read was actually "
    "observed splicing into it -- otherwise <code>te_exonized</code>. "
    "Neither is wrong. Do not read agreement between the two as "
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
        "id": "chimera_chimeric_reads_te_type",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Chimeric reads - TE type",
        "description": (
            "Per-sample chimeric junctions by TE type, from STAR's chimeric "
            "reads. The assembly screen uses these <strong>same words for "
            "a different measurement</strong> -- decided by transcript "
            "structure there, not junction direction; neither is wrong, but "
            "do not read agreement as corroboration (see Help)."
        ),
        "helptext": CROSS_SCREEN_NOTE.replace("<br><br>", ""),
    }
    if any(v for row in data.values() for v in row.values()):
        doc.update({
            "plot_type": "bargraph",
            "pconfig": {
                "id": "chimera_chimeric_reads_te_type_plot",
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
