#!/usr/bin/env python3
"""The report's guide to reading the gene-TE chimera evidence -- and nothing
more than a guide.

This pipeline does not rank chimera candidates.  It used to: a four-tier
confidence ladder, rendered here as the report's first chimera section.  The
ladder was removed because no experiment in this project established the
relative weight of its rungs, and the pipeline's own cohort analysis
contradicted the top one -- chimera_evidence_heatmap.py measured cross-screen
agreement sitting near its chance rate, the same class of result that had
already removed TElocal expression from the ladder (0d04e43).

Ranking on an unvalidated weighting is worse than not ranking, because a
tier column in a report is read as a verdict.  So the report now states what
is known about each line of evidence and stops.  Deciding which candidates
are real is a manual call made against the full catalogue,
results/chimera/candidates.tsv.gz.

Deliberately absent from this section: any candidate table.  Rendering even
an unranked top-N re-creates the thing that was removed, because whatever
order it happens to be in reads as importance.  Every pair is in the TSV,
with every evidence column; the reader sorts it for the question they have.
Guard 37 asserts no chimera section renders a table.

Emits two MultiQC custom-content documents:

  chimera_evidence_guide_mqc.json        html   how to weigh each signal
  chimera_evidence_composition_mqc.json  bar    pairs carrying each signal
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

PARENT_ID = "chimera"
PARENT_NAME = "Chimera"

# Flags emitted by chimera_evidence.py, in the order chimera_evidence.py
# builds them. Presentational only -- the whole point of this section is that
# no order among these is established.
FLAGS = [
    ("canonical", "Splice motif"),
    ("multi_sample", "Replicate support"),
    ("both_screens", "Called by both screens"),
    ("assembly_strand_match", "Assembly strand match"),
    ("telocal_expressed", "TE locus expressed"),
]

# (label, source tool, what the signal is, what THIS pipeline has measured
#  about it, standing)
#
# The source column exists because the labels alone are ambiguous: "splice
# motif" and "strand match" could plausibly come from either screen, and a
# reader weighing two signals needs to know whether they are independent
# measurements or two views of the same tool's output.
SIGNALS = [
    (
        "Splice motif",
        "STAR (chimeric junctions)",
        "A recognised splice motif on at least one junction "
        "(<code>canonical</code>).",
        "The best artifact discriminator available here. Real introns are "
        "~100% canonical, while template-switching, ligation and PCR chimeras "
        "carry no motif at all. A low overall rate is normal &mdash; what "
        "matters is enrichment within a donor group, not the absolute number.",
        "strong",
    ),
    (
        "Replicate support",
        "STAR (chimeric junctions)",
        "Seen in more than one sample (<code>multi_sample</code>).",
        "Weaker than it looks. A sequence-driven template switch recurs across "
        "libraries too, so recurrence does not separate a real chimera from a "
        "reproducible artifact.",
        "mixed",
    ),
    (
        "Called by both screens",
        "STAR + StringTie",
        "Found by the read-evidence <em>and</em> transcript-evidence screens "
        "(<code>both_screens</code>).",
        "Should be the strongest signal here &mdash; the two screens have "
        "opposite blind spots. Measured across a cohort it was not: agreement "
        "came out near its <strong>chance rate</strong>. Treat it as "
        "unresolved, and check the <strong>Evidence structure</strong> "
        "sections below for your own data before relying on it.",
        "unresolved",
    ),
    (
        "Assembly strand match",
        "StringTie (assembly)",
        "The assembled transcript's strand agrees with the gene's "
        "(<code>assembly_strand_match</code>).",
        "A consistency check on the assembly call rather than independent "
        "support. For <code>te_initiated</code> calls a mismatch usually means "
        "the gene hit is a spurious overlap, not real transcript connectivity.",
        "mixed",
    ),
    (
        "Read depth",
        "STAR (chimeric junctions)",
        "<strong>Not an evidence flag.</strong> Reported as "
        "<code>junction_reads</code> / <code>junction_events</code>.",
        "The metric most inflated by artifacts &mdash; a hot PCR chimera is "
        "often the deepest event in a run. Depth never promotes a pair here, "
        "and a high-depth row with no flags means exactly that.",
        "not-evidence",
    ),
    (
        "TE locus expressed",
        "TElocal",
        "Counted as an evidence flag. Reported as "
        "<code>telocal_count</code> when TElocal ran.",
        "One small 4-sample mouse experiment: 91% of junction-side pairs had "
        "an expressed locus, and the canonical rate was <em>lower</em> where "
        "it was (6.7% vs 10.2%, n&nbsp;=&nbsp;19,503). Too early to conclude "
        "anything from a single run &mdash; the correlation between "
        "junction-side pairs and TE locus expression needs testing properly. "
        "It stays a flag until that test exists &mdash; one small experiment "
        "is not enough to demote a signal.",
        "unresolved",
    ),
]

WEIGHT_STYLE = {
    "strong": ("#1a7f5a", "Best discriminator"),
    "mixed": ("#8a6d15", "Read with care"),
    "unresolved": ("#8a3a54", "Unresolved"),
    "not-evidence": ("#777", "Not evidence"),
}


def load(path):
    """Rows of a TSV as dicts, keyed by header name."""
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def guide_html(n_pairs, composition, n_no_flags):
    rows = []
    for label, source, what, measured, weight in SIGNALS:
        colour, badge = WEIGHT_STYLE[weight]
        rows.append(
            "<tr>"
            f'<td style="white-space:nowrap;"><strong>{label}</strong></td>'
            f'<td style="white-space:nowrap;color:#555;">{source}</td>'
            f'<td style="white-space:nowrap;color:{colour};">{badge}</td>'
            f"<td>{what}</td>"
            f"<td>{measured}</td>"
            "</tr>"
        )

    counted = ", ".join(
        f"<code>{flag}</code> {composition.get(flag, 0):,}" for flag, _ in FLAGS
    )
    return f"""
<p>What each column of the <strong>Candidates</strong> table above is worth,
and what this project has actually measured about it. Sort that table on the
signal your question needs &mdash; this is the reference for choosing which
one, and for knowing how far to trust it.</p>

<table class="table" style="width:100%; font-size: 90%;">
<thead><tr>
<th>Signal</th><th>Source</th><th>Standing</th><th>What it is</th>
<th>What has been measured about it</th>
</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>

<p style="margin-top:1em;">This run produced <strong>{n_pairs:,}</strong>
gene-TE pairs: {counted}. <strong>{n_no_flags:,}</strong> carry no evidence
flag at all.</p>

<p style="font-size: 85%; color: #888;">Full catalogue, one row per pair with
every evidence column:
<code>results/chimera/candidates.tsv.gz</code>.</p>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True,
                    help="chimera candidates.tsv.gz (chimera_evidence.py)")
    ap.add_argument("--out-guide", required=True)
    ap.add_argument("--out-composition", required=True)
    args = ap.parse_args()

    rows = list(load(args.evidence))
    n_pairs = len(rows)

    composition = {flag: 0 for flag, _ in FLAGS}
    n_no_flags = 0
    for r in rows:
        present = [f for f in r.get("evidence", ".").split(",") if f != "."]
        if not present:
            n_no_flags += 1
        for flag in present:
            if flag in composition:
                composition[flag] += 1

    guide_doc = {
        # NOT "chimera_evidence_guide": a section id must never equal a
        # parent_id. report_section_order's first pass matches MODULE anchors,
        # and a custom-content group's anchor is its parent_id -- so a section
        # sharing that name gets picked up by the module pass, whose order
        # semantics are inverted, and the whole group renders LAST. Measured.
        #
        # Every view in this group must still use the same parent_id, or
        # MultiQC silently splits the group into two report sections.
        "id": "chimera_signal_guide",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "section_name": "How to weigh this evidence",
        "description": (
            "What each line of chimera evidence is worth, and what this "
            "pipeline has actually measured about it. No ranking is produced."
        ),
        "plot_type": "html",
        "data": guide_html(n_pairs, composition, n_no_flags),
    }

    # Composition, not a ranking: how much of each signal the cohort produced.
    # Shown because "148 of 2,431 carry a splice motif" changes how the whole
    # table should be read, and is invisible from any individual row.
    # A bargraph whose every value is zero does not render as an empty plot --
    # MultiQC raises ValueError("No datasets to plot"), exits non-zero, and
    # NO REPORT IS WRITTEN AT ALL. A run that finds no gene-TE pairs is a
    # normal outcome (a clean library, or a species with a sparse TE
    # annotation), so it must not cost the user their entire report. Fall back
    # to an HTML section, the same way chimera_assembly_summary_mqc.py does
    # for its own empty case.
    counts_by_label = {
        **{label: composition[flag] for flag, label in FLAGS},
        "No evidence flag": n_no_flags,
    }
    if any(counts_by_label.values()):
        # One bar PER EVIDENCE TYPE, not one stacked bar split by type.
        # MultiQC stacks by default (stacking="relative"), which drew these as
        # segments of a single bar -- and that reads as a partition: mutually
        # exclusive slices summing to the cohort. It is the opposite of true.
        # A pair can carry all five flags at once, so the counts OVERLAP and
        # do not sum to anything meaningful. Giving each its own bar removes
        # the implied exclusivity.
        composition_body = {
            "plot_type": "bargraph",
            "pconfig": {
                "id": "chimera_evidence_composition_plot",
                "title": "Gene-TE pairs carrying each line of evidence",
                "ylab": "Gene-TE pairs",
                "cpswitch": False,
                "stacking": "group",
                # whole pairs: no decimals. tt_decimals is the key MultiQC
                # honours here; "format" is silently dropped.
                "tt_decimals": 0,
            },
            "categories": ["Gene-TE pairs"],
            "data": {
                label: {"Gene-TE pairs": count}
                for label, count in counts_by_label.items()
            },
        }
    else:
        composition_body = {
            "plot_type": "html",
            "data": (
                "<p>No gene-TE pairs were found in this run, so there is "
                "nothing to plot here. This is not an error: it means neither "
                "screen called a gene-TE chimera. The screens are independent "
                "of each other and of TE quantification, so TEcount and "
                "TElocal results are unaffected.</p>"
            ),
        }

    composition_doc = {
        "id": "chimera_evidence_composition",
        "parent_id": PARENT_ID,
        "parent_name": PARENT_NAME,
        "section_name": "Evidence composition",
        "description": (
            "How many gene-TE pairs carry each line of evidence. "
            "<strong>These bars overlap and do not sum to the cohort.</strong> "
            "A single pair can carry all five flags at once, so it is counted "
            "in several bars &mdash; they are independent counts, not slices "
            "of a whole, which is why they are drawn separately rather than "
            "stacked. Sources: splice motif and replicate support from STAR "
            "chimeric junctions, assembly strand match from StringTie, "
            "both-screens from the two together."
        ),
        **composition_body,
    }

    for path, doc in ((args.out_guide, guide_doc),
                      (args.out_composition, composition_doc)):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open_write(path) as fh:
            json.dump(doc, fh, indent=2)
            fh.write("\n")

    summary = ", ".join(f"{flag}: {composition[flag]}" for flag, _ in FLAGS)
    print(f"chimera evidence guide: {n_pairs} gene-TE pairs ({summary}, "
          f"no flag: {n_no_flags}) -> {args.out_guide}, {args.out_composition}")


if __name__ == "__main__":
    main()
