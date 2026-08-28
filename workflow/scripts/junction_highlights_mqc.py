#!/usr/bin/env python3
"""A "what to look at first" guide + ranked candidate table for the
chimera-JUNCTION screen, as MultiQC custom content.

The mirror of chimera_assembly_summary_mqc.py's highlights section, which
until now was the only guided entry point in the report -- and it belongs to
the assembly screen. Both screens are on by default, but the junction screen
is the one that always runs first and produces the larger table (tens of
thousands of events on real data), so a report with no ranking for it left
users staring at all_events.tsv.gz with no idea which rows matter.

Reads the merged gene<->TE table (chimera_counts.py's --out-te-events) rather
than the per-sample QC metrics, because ranking needs the per-event
annotation columns: canonical motif, replicate support, strand agreement and
(when telocal is enabled) whether the TE locus is actually expressed.

Ranking, strongest signal first:
  1. canonical == yes      -- a recognised splice motif.  Real introns are
                              ~100% canonical while template-switching,
                              ligation and PCR chimeras carry no motif, so
                              this is the single best artifact discriminator
                              in the screen.
  2. n_samples             -- reproducibility across the cohort.
  3. gene_strand_match     -- consistent transcript connectivity.
  4. total_reads           -- depth, deliberately last: read count alone is
                              the metric most inflated by artifacts (a hot
                              PCR chimera can be the deepest event in a run).

telocal_active is displayed but NOT ranked on.  It was originally the second
key, on the assumption that an independently expressed TE locus corroborates
the chimera.  On a real 4-sample mouse run it does the opposite: 91% of
gene-TE junctions had an expressed locus, and the canonical rate was LOWER
where the TE was expressed (6.7% at telocal_count > 10 vs 10.2% at <= 10,
n = 19,503).  A highly expressed locus yields more reads and so more chances
for template switching, so expression is context, not support.
"""
import argparse
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


def esc(value):
    return html.escape(str(value))


def load(path):
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--te-events", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=25)
    args = ap.parse_args()

    rows = list(load(args.te_events))
    n_total = len(rows)
    # telocal columns only exist when the telocal-annotation step ran; "." is
    # what chimera_counts writes for an absent column, so this also stays
    # correct if telocal is on but a run produced no annotated events.
    has_telocal = any(r.get("telocal_active", ".") in ("yes", "no") for r in rows)

    def rank_key(r):
        return (
            1 if r.get("canonical") == "yes" else 0,
            _int(r.get("n_samples")),
            1 if r.get("gene_strand_match") == "yes" else 0,
            _int(r.get("total_reads")),
        )

    top = sorted(rows, key=rank_key, reverse=True)[: args.top_n]

    n_canonical = sum(1 for r in rows if r.get("canonical") == "yes")
    n_active = sum(1 for r in rows if r.get("telocal_active") == "yes")
    n_multi = sum(1 for r in rows if _int(r.get("n_samples")) > 1)

    cols = ["event_id", "chimera_type", "gene_id", "te_id", "te_subfamily",
            "te_class", "canonical", "gene_strand_match"]
    if has_telocal:
        cols.append("telocal_active")
    cols += ["n_samples", "total_reads"]
    col_label = {
        "event_id": "Event", "chimera_type": "Class", "gene_id": "Gene",
        "te_id": "TE insertion", "te_subfamily": "Subfamily",
        "te_class": "TE class", "canonical": "Splice motif",
        "gene_strand_match": "Strand match", "telocal_active": "TE expressed",
        "n_samples": "Samples", "total_reads": "Reads",
    }

    thead = "".join(f"<th>{esc(col_label[c])}</th>" for c in cols)
    trows = []
    for r in top:
        cells = []
        for c in cols:
            value = r.get(c, ".")
            if c == "canonical" and value == "yes":
                cells.append(f"<td>{esc(value)} &#9733;</td>")
            else:
                cells.append(f"<td>{esc(value)}</td>")
        trows.append("<tr>" + "".join(cells) + "</tr>")
    if not trows:
        trows = [f"<tr><td colspan={len(cols)}>No gene-TE junctions.</td></tr>"]

    pct = (100.0 * n_canonical / n_total) if n_total else 0.0
    pct_active = (100.0 * n_active / n_total) if n_total else 0.0
    telocal_line = (
        f"<li><strong>{n_active}</strong> ({pct_active:.0f}%) involve a TE "
        "locus TElocal finds expressed (<code>telocal_active</code>). Read "
        "this as <em>context, not support</em>: on real data an expressed "
        "locus is the common case and its junctions are <em>less</em> likely "
        "to carry a splice motif, since more reads means more chances for "
        "template switching. It is shown, not ranked on.</li>"
        if has_telocal else
        "<li>TElocal is disabled, so no per-locus expression context is "
        "available for the TE side.</li>"
    )

    body = f"""
<p>This screen reads STAR's chimeric junctions: reads that cannot be
explained by one linear alignment. It is annotation-blind by construction,
so it catches breakpoints no assembly would predict -- but it also cannot
see a gene-TE chimera spliced through an ordinary canonical intron, which
is exactly what the <em>assembly</em> screen is for. Neither screen alone
is complete; an event found by both is the strongest evidence available
here.</p>

<p><strong>What to look at first:</strong></p>
<ul>
<li><strong>{n_canonical} of {n_total}</strong> gene-TE junctions
({pct:.1f}%) carry a recognised splice motif
(<code>canonical: yes</code>, marked &#9733; below). Start here: real introns
are ~100% canonical, while template-switching, ligation and PCR chimeras
carry no motif at all. A low overall percentage is normal -- what matters is
enrichment, and the honest comparison is <em>within a donor group</em>:
<code>gene_to_te</code> against <code>gene_to_gene</code>/
<code>gene_to_other</code>, not against <code>other</code>. See the
"Canonical rate by direction" plot.</li>
{telocal_line}
<li><strong>{n_multi}</strong> events are seen in more than one sample.
Prefer reproducibility over depth -- though note it is weaker than a splice
motif here, because a sequence-driven template switch recurs across
libraries too.</li>
<li><code>direction_ambiguous: yes</code> means the breakpoint hit a gene
<em>and</em> a TE (common -- TEs sit inside gene bodies), so the reported
direction was decided by branch order, not by evidence. Treat the call as
provisional.</li>
<li><code>te_id</code> is the individual insertion and joins against TElocal
rows; <code>te_subfamily</code> joins against TEcount rows. If a subfamily
looks expressed in TEcount, this is how you find which copy.</li>
</ul>

<table class="table" style="width:100%; font-size: 90%;">
<thead><tr>{thead}</tr></thead>
<tbody>
{"".join(trows)}
</tbody>
</table>
<p style="font-size: 85%; color: #888;">Showing the top {len(top)} of
{n_total} gene-TE junctions (ranked: splice motif first, then replicate
support, then strand match, then reads). Full table:
<code>results/chimera_junction/te-gene-chimeras.tsv.gz</code>.</p>
"""

    doc = {
        # parent_id must match junction_qc_mqc.py's "chimera" group exactly
        # -- a mismatch splits the Chimera section in two in the report.
        "id": "junction_highlights",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "What to look at",
        "description": (
            "How to read the chimera-junction output, and this run's "
            "highest-confidence gene-TE junctions."
        ),
        "plot_type": "html",
        "data": body,
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(
        f"junction highlights: {n_total} gene-TE junctions, "
        f"{n_canonical} canonical -> {args.out}"
    )


if __name__ == "__main__":
    main()
