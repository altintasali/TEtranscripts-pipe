#!/usr/bin/env python3
"""Build MultiQC custom-content for the transcript-evidence (assembly) screen:
  chimera_assembly_classes_mqc.json  candidate counts by chimera_type,
                                      split by cross-confirmation with the
                                      read-evidence screen when that column
                                      is present -- otherwise just totals
                                      per class.
  chimera_assembly_highlights_mqc.json  what this screen structurally can and
                                      cannot see, plus the counts that
                                      qualify its own output.

This section used to also rank and render this run's top candidates, keyed
(junction-confirmed, strand-matched, highest TPM).  It was one of three
rankings shipping in the same report, and it led on the one signal
chimera_evidence_heatmap.py measured at roughly its chance rate.  All three
are gone -- the pipeline no longer ranks candidates at all (see
chimera_evidence_guide_mqc.py); dropping the table here also dropped this
script's only use of tpm_matrix.tsv.gz, so it is no longer an input.

Reads the candidates table (or, when present, the cross-referenced
transcripts_with_read_support.tsv.gz -- pass whichever is available as
--candidates).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gz_io import open_read, open_write

CLASS_ORDER = [
    "te_initiated", "te_initiated_intergenic", "te_exonized",
    "te_terminated", "unspliced_te_only",
]

CLASS_LABEL = {
    "te_initiated": "TE-initiated",
    "te_initiated_intergenic": "TE-initiated (no gene match)",
    "te_exonized": "TE-exonized",
    "te_terminated": "TE-terminated",
    "unspliced_te_only": "Unspliced (low confidence)",
}


def read_table(path):
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True,
                     help="candidates.tsv.gz, or candidates_with_junction_evidence.tsv.gz if available")
    ap.add_argument("--out-classes", required=True)
    ap.add_argument("--out-highlights", required=True)
    ap.add_argument("--out-strand-rate", required=True)
    args = ap.parse_args()

    header, rows = read_table(args.candidates)
    has_cross_evidence = "confirmed_by_junction_screen" in header

    # --- classes barplot -----------------------------------------------
    if has_cross_evidence:
        counts = {c: {"confirmed": 0, "unconfirmed": 0} for c in CLASS_ORDER}
        for r in rows:
            c = r.get("chimera_type", "")
            if c not in counts:
                continue
            key = "confirmed" if r.get("confirmed_by_junction_screen") == "yes" else "unconfirmed"
            counts[c][key] += 1
        data_labels_note = (
            "confirmed = also found by the chimera-junction (read-level) screen -- "
            "independent evidence, higher confidence."
        )
    else:
        counts = {c: {"count": 0} for c in CLASS_ORDER}
        for r in rows:
            c = r.get("chimera_type", "")
            if c in counts:
                counts[c]["count"] += 1
        data_labels_note = (
            "chimera.junction is disabled, so these counts have no independent "
            "cross-confirmation -- enable it too for higher-confidence calls."
        )

    classes_doc = {
        "id": "chimera_assembly_classes",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Assembly - composition by class",
    }
    # counts values are {"count": N} or {"confirmed": N, "unconfirmed": N}
    # depending on whether the junction screen ran, so sum values rather
    # than assuming either key.
    if any(sum(v.values()) for v in counts.values()):
        classes_doc.update({
            "description": (
                "StringTie-assembly chimera candidates by chimera_type. "
                + data_labels_note
            ),
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_assembly_classes_plot",
                "title": "Chimera-assembly candidates by class",
                "ylab": "candidates",
                "cpswitch": False,
                "use_legend": True,
            },
            "data": {CLASS_LABEL.get(c, c): v for c, v in counts.items()},
        })
    else:
        # MultiQC's bargraph raises "No datasets to plot" on an all-zero
        # bar document and takes the whole report down with it, so document
        # the empty result as html instead of shipping an unplottable bar.
        # A run with no assembled TE candidates is normal (small genome, low
        # depth, or genuinely nothing there) -- it must not fail the report.
        classes_doc.update({
            "description": "StringTie-assembly chimera candidates by chimera_type.",
            "plot_type": "html",
            "data": (
                "<p>No gene-TE chimera candidates were assembled in this run, "
                "so there is nothing to plot here. This is not an error: it "
                "means StringTie found no transcript whose structure implicates "
                "a TE. The junction screen is independent and may still have "
                "calls.</p>"
            ),
        })
    os.makedirs(os.path.dirname(args.out_classes) or ".", exist_ok=True)
    with open_write(args.out_classes) as fh:
        json.dump(classes_doc, fh, indent=2)
        fh.write("\n")

    # --- strand-match rate by class -------------------------------------
    # The structural mirror of the read screen's splice-motif rate: a per-class
    # quality rate rather than a per-class count. Derived from columns already
    # in the candidates table, so it costs nothing extra.
    #
    # Deliberately NOT mirrored: the read screen's "gene-TE subset". That shows
    # what share of all chimeric junctions involve a TE; classify_chimera_
    # assembly.py only ever emits TE-involving classes, so the assembly
    # equivalent would be a constant 100% and would say nothing.
    strand_rate, strand_counts = {}, {}
    for cls in CLASS_ORDER:
        in_cls = [r for r in rows if r.get("chimera_type") == cls]
        matched = sum(1 for r in in_cls if r.get("strand_match") == "yes")
        strand_rate[cls] = round(100.0 * matched / len(in_cls), 1) if in_cls else 0.0
        strand_counts[cls] = matched

    if any(strand_counts.values()):
        # One bar PER CLASS. The first shape here put every class as a
        # category of a single "strand match" bar, and MultiQC stacks by
        # default -- so the percentages were summed and the axis ran past
        # 100%. Percentages of different denominators must never share a
        # stack; each class is its own bar with its own 0-100 scale.
        strand_body = {
            "plot_type": "bar",
            "pconfig": {
                "id": "chimera_assembly_strand_rate_plot",
                "title": "Strand-match rate by chimera class",
                "ylab": "% of candidates in the class",
                # Raw counts first, then the rate. ymax belongs to the
                # PERCENTAGE dataset only -- at plot level it also capped the
                # counts view, clipping any class with more than 100
                # candidates. cpswitch stays off: a strand-match RATE is
                # matched/total within a class, not a share of the whole, so
                # MultiQC's own percentage would be a different number.
                "cpswitch": False,
                "stacking": "group",
                "data_labels": [
                    {"name": "Strand-matched candidates", "tt_decimals": 0},
                    {"name": "% strand-matched", "tt_decimals": 1, "ymax": 100},
                ],
            },
            "data": [
                {cls: {"Strand-matched candidates": strand_counts[cls]}
                 for cls in CLASS_ORDER},
                {cls: {"% strand-matched": strand_rate[cls]} for cls in CLASS_ORDER},
            ],
        }
    else:
        # all-zero plots crash MultiQC outright -- see chimera_evidence_guide_mqc
        strand_body = {
            "plot_type": "html",
            "data": ("<p>No candidate in any class has a matching strand, so "
                     "there is nothing to plot. With few candidates this is "
                     "unremarkable; across a large assembly it suggests the "
                     "gene hits are spurious overlaps.</p>"),
        }

    strand_doc = {
        "id": "chimera_assembly_strand_rate",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Assembly - strand-match rate by class",
        "description": (
            "<em>What \"strand match\" means:</em> StringTie assembles each "
            "transcript on a strand, and the gene it overlaps is annotated on "
            "a strand. <strong>Strand match = those two agree.</strong> They "
            "have to, for the transcript to actually be a fusion of that gene "
            "with that TE: a transcript on the opposite strand from the gene "
            "is running the other way and cannot be transcribed from it, so "
            "the overlap is coincidental &mdash; the two features merely sit "
            "in the same place in the genome. "
            "<br><br><em>How to read it:</em> each bar is one chimera class, "
            "scored independently, so the bars do not sum to anything. A "
            "class with a low rate is one to distrust: for "
            "<code>te_initiated</code> especially, a mismatch usually means "
            "the gene hit is a spurious overlap rather than real transcript "
            "connectivity. This is a consistency check on the assembly, not "
            "independent support &mdash; the read screen\'s splice-motif "
            "rate is the closer thing to evidence."
        ),
        **strand_body,
    }
    os.makedirs(os.path.dirname(args.out_strand_rate) or ".", exist_ok=True)
    with open_write(args.out_strand_rate) as fh:
        json.dump(strand_doc, fh, indent=2)
        fh.write("\n")

    # --- screen notes: blind spot + qualifying counts --------------------
    # This section used to render its own ranked top-N, keyed
    # (junction-confirmed, strand-matched, TPM). That was one of three
    # rankings shipping in the same report, and it led on the very signal
    # chimera_evidence_heatmap.py measured at roughly its chance rate. The
    # pipeline no longer ranks candidates anywhere; what stays here is what
    # only this screen can say about itself.
    n_total = len(rows)
    n_confirmed = (
        sum(1 for r in rows if r.get("confirmed_by_junction_screen") == "yes")
        if has_cross_evidence else None
    )
    n_strand_ok = sum(1 for r in rows if r.get("strand_match") == "yes")
    n_unspliced = sum(
        1 for r in rows if r.get("chimera_type") == "unspliced_te_only"
    )

    confirmed_line = (
        f"<li><strong>{n_confirmed} of {n_total}</strong> candidates are also "
        "found by the read-evidence screen. The two methods have opposite "
        "blind spots, so agreement should be meaningful &mdash; but measured "
        "across a cohort it has come out near its <strong>chance rate</strong>. "
        "It carries no more weight than any other flag here; see "
        "the <strong>Evidence structure</strong> sections below for your "
        "own data.</li>"
        if has_cross_evidence else
        "<li>chimera.junction is currently disabled, so no cross-confirmation "
        "is available &mdash; enabling it adds an independent evidence source "
        "for these same candidates.</li>"
    )

    html = f"""
<p><strong>What this screen sees.</strong> Gene-TE chimeras inferred from
StringTie assembly structure: TE-initiated, exonized and TE-terminated
transcripts spliced through an ordinary, canonical intron.</p>

<p><strong>What it cannot see.</strong> Any structure an assembler would not
build &mdash; the breakpoints that only show up as reads STAR cannot explain
with one linear alignment. That gap is what the read-evidence screen covers,
which is why the two are kept separate rather than merged.</p>

<p><strong>Qualifying this screen's output:</strong></p>
<ul>
{confirmed_line}
<li><strong>{n_strand_ok} of {n_total}</strong> candidates have
<code>strand_match: yes</code>. For <code>te_initiated</code> calls a mismatch
usually means the "gene" hit is a spurious overlap, not real transcript
connectivity.</li>
<li><strong>{n_unspliced}</strong> are <code>unspliced_te_only</code>, which
carry no splice evidence at all &mdash; useful for manual triage, not for a
confident call on their own.</li>
<li>Check <code>te_hits_all</code> / <code>gene_hits_all</code> in the full
table: more than one entry means the locus is ambiguous (a multi-copy TE
family, or overlapping gene isoforms) &mdash; the single <code>te_id</code>/
<code>matched_gene_id</code> reported is just the first hit, not necessarily
the right one.</li>
<li>Candidates are structural calls from the merged assembly, not per-sample.
See <code>tpm_matrix.tsv.gz</code> for real expression, and require support
across multiple replicates of a condition before trusting a call as biological
rather than assembly noise.</li>
</ul>

<p style="font-size: 85%; color: #888;">This screen's calls are merged with
the read-evidence screen's into the <strong>Candidates</strong> table
above.</p>
"""

    highlights_doc = {
        "id": "chimera_assembly_highlights",
        "parent_id": "chimera",
        "parent_name": "Chimera",
        "section_name": "Assembly - what this screen sees",
        "description": (
            "The transcript-evidence screen's blind spot, and the counts that "
            "qualify its output."
        ),
        "plot_type": "html",
        "data": html,
    }
    os.makedirs(os.path.dirname(args.out_highlights) or ".", exist_ok=True)
    with open_write(args.out_highlights) as fh:
        json.dump(highlights_doc, fh, indent=2)
        fh.write("\n")

    print(
        f"chimera_assembly summary: {n_total} candidates"
        + (f", {n_confirmed} junction-confirmed" if has_cross_evidence else "")
        + f" -> {args.out_classes}, {args.out_highlights}"
    )


if __name__ == "__main__":
    main()
