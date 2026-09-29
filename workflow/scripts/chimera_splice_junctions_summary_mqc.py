#!/usr/bin/env python3
"""Build MultiQC custom-content for the SJ.out.tab junction screen
(chimera.splice_junctions):
  chimera_splice_junctions_highlights_mqc.json  per-sample gene<->TE junction
                                      counts by direction, plus cohort totals.
  chimera_splice_junctions_te_type_mqc.json     cohort-wide candidate counts
                                      by chimera_type, with te_initiated split
                                      into upstream/internal.

Modelled on chimera_assembly_summary_mqc.py: like the assembly screen, this
one has no per-candidate row worth rendering individually here (that is what
the shared Candidates table is for), so both sections are aggregate counts.

Before this script existed, this screen -- unlike chimeric_reads and
assembly -- had NO always-on report presence at all: its only two views
(splice_junctions_pca_*/heatmap_*) are gated behind
chimera.splice_junctions.qc.enabled, which defaults to false. A run with the
screen enabled but that flag left at its default showed nothing for it, even
though sj_canonical alone flagged far more pairs than cr_canonical in a
recent cohort (15,913 vs 1,647) and the "TE analysis" overview tells the
reader all three screens are running.

No new computation: both inputs are already-classified tables
(classify_chimera_splice_junctions.py / chimera_splice_junctions_counts.py)
-- this only regroups columns that already exist.

Reads:
  --te-events   results/chimera/splice_junctions/te-gene-junctions.tsv.gz
                (cohort-level: one row per UNIQUE gene<->TE junction, already
                deduplicated across samples -- chimera_splice_junctions_
                counts.py's --out-te-events). Used for the TE-type
                composition (section 2) and the cohort totals quoted in
                section 1's description.
  --per-sample  one results/chimera/splice_junctions/per_sample/{sample}_
                junctions_te-gene-junctions.tsv.gz per sample (same columns,
                NOT deduplicated across samples -- classify_chimera_splice_
                junctions.py's --te-out). Used for section 1's per-sample
                bar chart: summing the cohort-level file's rows per sample
                is not possible (it carries no sample column, by design --
                see chimera_splice_junctions_counts.py's module docstring),
                and summing occurrence counts from the per-sample files
                would double-count any junction seen in more than one
                sample if used for the cohort totals instead.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

DIRECTIONS = ["gene_to_te", "te_to_gene"]

# te_initiated is split by te_initiated_detail (classify_chimera_splice_
# junctions.py) into upstream/internal -- the distinction only this screen
# makes, and the point of the classification fix that added the column:
# see that script's module docstring for the full bug/fix reasoning.
# te_exonized/te_terminated have no such split (te_initiated_detail is "."
# for gene_to_te events).
TE_TYPE_CATEGORIES = [
    "te_initiated_upstream", "te_initiated_internal",
    "te_exonized", "te_terminated", "antisense_to_gene",
]
TE_TYPE_LABEL = {
    "te_initiated_upstream": "TE-initiated (upstream)",
    "te_initiated_internal": "TE-initiated (internal, skips an earlier exon)",
    "te_exonized": "TE-exonized",
    "te_terminated": "TE-terminated",
    "antisense_to_gene": "Antisense to the gene (opposite strand)",
}


def read_table(path):
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]
    return rows


def te_type_key(row):
    """chimera_type, with te_initiated split by te_initiated_detail."""
    c = row.get("chimera_type", "")
    if c == "te_initiated":
        detail = row.get("te_initiated_detail", ".")
        if detail in ("upstream", "internal"):
            return f"te_initiated_{detail}"
        return None  # defensive: classify_chimera_splice_junctions.py
        # always sets one of the two for a real te_initiated row.
    if c in ("te_exonized", "te_terminated", "antisense_to_gene"):
        return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--te-events", required=True,
                     help="results/chimera/splice_junctions/te-gene-junctions.tsv.gz")
    ap.add_argument("--per-sample", required=True, nargs="+",
                     help="results/chimera/splice_junctions/per_sample/"
                     "{sample}_junctions_te-gene-junctions.tsv.gz, one per sample")
    ap.add_argument("--sample-names", required=True, nargs="+")
    ap.add_argument("--out-highlights", required=True)
    ap.add_argument("--out-te-type", required=True)
    args = ap.parse_args()

    if len(args.per_sample) != len(args.sample_names):
        sys.exit("error: --per-sample and --sample-names must have equal length")

    cohort_rows = read_table(args.te_events)

    # --- section 1: per-sample direction composition + cohort totals ----
    per_sample_counts = {}
    for path, sample in zip(args.per_sample, args.sample_names):
        counts = dict.fromkeys(DIRECTIONS, 0)
        for r in read_table(path):
            d = r.get("direction", "")
            if d in counts:
                counts[d] += 1
        per_sample_counts[sample] = counts

    n_cohort_total = len(cohort_rows)
    n_cohort_by_dir = {
        d: sum(1 for r in cohort_rows if r.get("direction") == d)
        for d in DIRECTIONS
    }

    highlights_doc = {
        "id": "chimera_splice_junctions_highlights",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Splice junctions - what this screen sees",
    }
    if any(sum(c.values()) for c in per_sample_counts.values()):
        highlights_doc.update({
            "description": (
                f"Per-sample gene↔TE junctions from STAR's own splice "
                f"junctions (SJ.out.tab), by direction. This cohort "
                f"produced <strong>{n_cohort_total:,}</strong> unique "
                f"gene↔TE junctions ({n_cohort_by_dir['gene_to_te']:,} "
                f"gene_to_te, {n_cohort_by_dir['te_to_gene']:,} te_to_gene). "
                f"Canonical rate isn't plotted here: it's 100% by "
                f"construction under the default config "
                f"(require_canonical: true admits only canonical "
                f"junctions) -- if you've set require_canonical: false, "
                f"check the canonical column in te-gene-junctions.tsv.gz "
                f"directly instead."
            ),
            "helptext": (
                "<p><strong>What this screen sees.</strong> A TE splicing "
                "into a gene (or a gene splicing into a TE) through an "
                "ordinary, canonical intron, straight from STAR's normal "
                "splice-junction output (SJ.out.tab) -- no extra STAR pass "
                "and no assembly step, unlike the other two screens.</p>"
                "<p><strong>What it cannot see.</strong> Anything requiring "
                "a chimeric (non-linear) alignment -- that is the "
                "chimeric-reads screen's job -- or a fully assembled "
                "multi-exon transcript structure -- that is the assembly "
                "screen's job.</p>"
                "<p>Counts here are per sample and NOT deduplicated across "
                "samples (the same junction seen in two samples counts "
                "twice here); the cohort totals quoted above ARE "
                "deduplicated, one row per unique junction "
                "(te-gene-junctions.tsv.gz).</p>"
                "<p style=\"font-size: 85%; color: #888;\">This screen's "
                "calls are merged with the other two screens' into the "
                "<strong>Candidates</strong> table above.</p>"
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_splice_junctions_highlights_plot",
                "title": "Splice-junction gene-TE events by direction",
                "ylab": "junctions",
                "cpswitch": True,
                "cpswitch_counts_label": "Junction counts",
                "cpswitch_percent_label": "% of this sample's gene-TE junctions",
                "tt_decimals": 0,
            },
            "data": per_sample_counts,
        })
    else:
        # an all-zero bargraph makes MultiQC exit non-zero and write no
        # report -- see chimera_evidence_guide_mqc.py for the same trap.
        highlights_doc.update({
            "description": (
                "Per-sample gene↔TE junctions from STAR's own splice "
                "junctions (SJ.out.tab)."
            ),
            "plot_type": "html",
            "data": (
                "<p>No gene↔TE splice junctions were found in this "
                "run, so there is nothing to plot here. This is not an "
                "error: it means no junction in any sample's SJ.out.tab "
                "had one breakpoint in a gene and the other in a TE. The "
                "other two screens are independent and may still have "
                "calls.</p>"
            ),
        })
    os.makedirs(os.path.dirname(args.out_highlights) or ".", exist_ok=True)
    with open_write(args.out_highlights) as fh:
        json.dump(highlights_doc, fh, indent=2)
        fh.write("\n")

    # --- section 2: cohort-wide TE-type composition ----------------------
    te_type_counts = dict.fromkeys(TE_TYPE_CATEGORIES, 0)
    for r in cohort_rows:
        key = te_type_key(r)
        if key in te_type_counts:
            te_type_counts[key] += 1

    te_type_doc = {
        "id": "chimera_splice_junctions_te_type",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Splice junctions - TE type",
    }
    if any(te_type_counts.values()):
        te_type_doc.update({
            "description": (
                "Cohort-wide gene-TE junctions by chimera_type. "
                "<code>te_initiated</code> is split into "
                "<em>upstream</em>/<em>internal</em> by "
                "<code>te_initiated_detail</code> -- a distinction only "
                "this screen makes (see Help)."
            ),
            "helptext": (
                "<code>te_initiated_upstream</code>: the acceptor exon is "
                "the gene's own most-5' annotated exon -- a conventional "
                "TE-initiated transcript. <code>te_initiated_internal</code>"
                ": the gene has an annotated exon further upstream that "
                "this junction's transcript skips -- an intronic TE "
                "promoter. Both used to be reported as a single "
                "<code>te_initiated</code> or <code>te_exonized</code> "
                "call depending on where the TE's coordinates happened to "
                "sit relative to the gene's overall span, which is a bug "
                "fixed in classify_chimera_splice_junctions.py: an "
                "intronic TE promoter splicing past the gene's real first "
                "exon (now <code>te_initiated_internal</code>) used to be "
                "scored <code>te_exonized</code> just because its "
                "coordinates fell inside the gene body, even though the "
                "junction itself says the transcript starts there -- see "
                "that script's module docstring for the full reasoning. "
                "<code>te_exonized</code> and <code>te_terminated</code> "
                "keep their usual meaning: an internal exon, or the last "
                "one, respectively."
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_splice_junctions_te_type_plot",
                "title": "Splice-junction gene-TE events by TE type",
                "ylab": "junctions",
                # One category per bar (a single cohort-wide count), so a
                # percentage view would read 100% for every class -- same
                # reasoning as chimera_assembly_summary_mqc.py's classes_doc.
                "cpswitch": False,
                "use_legend": True,
                # BUG FIXED 2026: without this, every bar showed "1.00" /
                # "0.00" instead of "1" / "0". MultiQC computes this plot's
                # hoverformat from whether the category values are int or
                # float (bargraph.py); the JSON round-trip through custom
                # content silently promotes plain ints to floats, so the
                # auto-detect always picked ",.2f" here. chimera_chimeric_
                # reads_te_type_mqc.py's per-sample plot already sets this
                # for the same reason -- these are whole junction counts,
                # never fractional.
                "tt_decimals": 0,
            },
            "data": {TE_TYPE_LABEL[c]: {"count": n}
                     for c, n in te_type_counts.items()},
        })
    else:
        te_type_doc.update({
            "description": "Cohort-wide gene-TE junctions by chimera_type.",
            "plot_type": "html",
            "data": (
                "<p>No gene↔TE splice junctions were found in this "
                "run, so there is nothing to plot here. This is not an "
                "error: it means no junction in any sample's SJ.out.tab "
                "had one breakpoint in a gene and the other in a TE.</p>"
            ),
        })
    os.makedirs(os.path.dirname(args.out_te_type) or ".", exist_ok=True)
    with open_write(args.out_te_type) as fh:
        json.dump(te_type_doc, fh, indent=2)
        fh.write("\n")

    print(
        f"chimera_splice_junctions summary: {n_cohort_total} unique "
        f"gene-TE junctions across {len(args.sample_names)} sample(s) "
        f"-> {args.out_highlights}, {args.out_te_type}"
    )


if __name__ == "__main__":
    main()
