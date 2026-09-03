#!/usr/bin/env python3
"""What the read-evidence (chimeric-junction) screen can and cannot see, plus
the screen-level counts that qualify its own output.

This section used to also render a ranked top-25, keyed on
(canonical, n_samples, gene_strand_match, total_reads).  That was one of three
rankings shipping in the same report -- the confidence ladder in
chimera_evidence.py sorted differently, and the assembly screen's highlights
section sorted differently again -- so a reader's candidate list depended on
which section they happened to scroll to.  All three are gone: the pipeline no
longer ranks candidates at all (see chimera_evidence_guide_mqc.py), and this
section keeps only what is genuinely specific to this screen: its blind spot,
and the counts a reader needs to judge its output.

Reads the merged gene<->TE table (chimera_reads_counts.py's --out-te-events) because
the qualifying counts are per-event annotation columns: canonical motif,
replicate support, and (when telocal is enabled) whether the TE locus is
expressed.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write


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
    args = ap.parse_args()

    rows = list(load(args.te_events))
    n_total = len(rows)
    # telocal columns only exist when the telocal-annotation step ran; "." is
    # what chimera_reads_counts writes for an absent column, so this also stays
    # correct if telocal is on but a run produced no annotated events.
    has_telocal = any(r.get("telocal_active", ".") in ("yes", "no") for r in rows)

    n_canonical = sum(1 for r in rows if r.get("canonical") == "yes")
    n_active = sum(1 for r in rows if r.get("telocal_active") == "yes")
    n_multi = sum(1 for r in rows if _int(r.get("n_samples")) > 1)

    pct = (100.0 * n_canonical / n_total) if n_total else 0.0
    pct_active = (100.0 * n_active / n_total) if n_total else 0.0

    telocal_line = (
        f"<li><strong>{n_active}</strong> ({pct_active:.0f}%) involve a TE "
        "locus TElocal finds expressed (<code>telocal_active</code>). Read "
        "this as <em>context, not support</em>: on real data an expressed "
        "locus is the common case and its junctions are <em>less</em> likely "
        "to carry a splice motif, since more reads means more chances for "
        "template switching.</li>"
        if has_telocal else
        "<li>TElocal is disabled, so no per-locus expression context is "
        "available for the TE side.</li>"
    )

    body = f"""
<p><strong>What this screen sees.</strong> STAR's chimeric junctions: reads
that cannot be explained by one linear alignment. It is annotation-blind by
construction, so it catches breakpoints no assembler would predict.</p>

<p><strong>What it cannot see.</strong> A gene-TE chimera spliced through an
ordinary, canonical intron &mdash; that read aligns linearly, so it never
reaches this screen at all. That gap is exactly what the transcript-evidence
(assembly) screen covers, which is why the two are kept separate rather than
merged.</p>

<p><strong>Qualifying this screen's output:</strong></p>
<ul>
<li><strong>{n_canonical} of {n_total}</strong> gene-TE junctions
({pct:.1f}%) carry a recognised splice motif (<code>canonical: yes</code>).
A low overall percentage is normal &mdash; what matters is enrichment, and the
honest comparison is <em>within a donor group</em>: <code>gene_to_te</code>
against <code>gene_to_gene</code>/<code>gene_to_other</code>, not against
<code>other</code>. See "splice-motif rate by junction class" below.</li>
{telocal_line}
<li><strong>{n_multi}</strong> events are seen in more than one sample. Note
this is weaker than a splice motif: a sequence-driven template switch recurs
across libraries too.</li>
<li><code>direction_ambiguous: yes</code> means the breakpoint hit a gene
<em>and</em> a TE (common &mdash; TEs sit inside gene bodies), so the reported
direction was decided by branch order, not by evidence. Treat the call as
provisional.</li>
<li><code>te_id</code> is the individual insertion and joins against TElocal
rows; <code>te_subfamily</code> joins against TEcount rows. If a subfamily
looks expressed in TEcount, this is how you find which copy.</li>
</ul>

<p style="font-size: 85%; color: #888;">This screen's calls are merged with
the assembly screen's into the <strong>Candidates</strong> table above.</p>
"""

    doc = {
        # parent_id must match chimera_reads_qc_mqc.py and sample_qc.R's chimera
        # view EXACTLY -- every view of this screen shares one group, and a
        # mismatch silently splits it into two report sections.
        "id": "chimera_reads_highlights",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Reads - what this screen sees",
        "description": (
            "The read-evidence screen's blind spot, and the counts that "
            "qualify its output."
        ),
        "plot_type": "html",
        "data": body,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open_write(args.out) as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(
        f"junction screen notes: {n_total} gene-TE junctions, "
        f"{n_canonical} canonical -> {args.out}"
    )


if __name__ == "__main__":
    main()
