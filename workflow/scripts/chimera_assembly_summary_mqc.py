#!/usr/bin/env python3
"""Build MultiQC custom-content for the chimera-assembly screen:
  chimera_assembly_classes_mqc.json  candidate counts by chimera_type,
                                      split by cross-confirmation with the
                                      junction screen when that column is
                                      present (results/chimera_assembly/
                                      candidates_with_junction_evidence.tsv.gz)
                                      -- otherwise just totals per class.
  chimera_assembly_highlights_mqc.json  a written "what to look at first"
                                      guide plus a ranked table of this
                                      run's actual top candidates (highest
                                      confidence first: cross-confirmed,
                                      strand-matched, highest expression).

Reads results/chimera_assembly/candidates.tsv.gz (or, when present,
candidates_with_junction_evidence.tsv.gz -- pass whichever is available as
--candidates) and results/chimera_assembly/tpm_matrix.tsv.gz.
"""
import argparse
import gzip
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


def load_max_tpm(path):
    """{transcript_id: max TPM across samples}."""
    max_tpm = {}
    with open_read(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            tid = parts[0]
            vals = [float(v) for v in parts[1:] if v]
            max_tpm[tid] = max(vals) if vals else 0.0
    return max_tpm


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True,
                     help="candidates.tsv.gz, or candidates_with_junction_evidence.tsv.gz if available")
    ap.add_argument("--tpm-matrix", required=True)
    ap.add_argument("--out-classes", required=True)
    ap.add_argument("--out-highlights", required=True)
    ap.add_argument("--top-n", type=int, default=25)
    args = ap.parse_args()

    header, rows = read_table(args.candidates)
    has_cross_evidence = "confirmed_by_junction_screen" in header
    max_tpm = load_max_tpm(args.tpm_matrix)

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
        "parent_id": "chimera_transcripts",
        "parent_name": "Gene-TE chimeras: transcript evidence",
        "section_name": "Candidates by class",
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
    os.makedirs(os.path.dirname(args.out_classes), exist_ok=True)
    with open_write(args.out_classes) as fh:
        json.dump(classes_doc, fh, indent=2)
        fh.write("\n")

    # --- highlights: guidance + ranked top-N table ----------------------
    def rank_key(r):
        confirmed = 1 if (has_cross_evidence and r.get("confirmed_by_junction_screen") == "yes") else 0
        strand_ok = 1 if r.get("strand_match") == "yes" else 0
        tpm = max_tpm.get(r.get("transcript_id", ""), 0.0)
        return (confirmed, strand_ok, tpm)

    ranked = sorted(rows, key=rank_key, reverse=True)
    top = ranked[: args.top_n]

    cols = ["transcript_id", "chimera_type", "te_id", "te_family", "te_class",
            "matched_gene_id", "strand_match"]
    if has_cross_evidence:
        cols.append("confirmed_by_junction_screen")
    cols.append("max_tpm")

    col_label = {
        "transcript_id": "Transcript", "chimera_type": "Class", "te_id": "TE",
        "te_family": "Family", "te_class": "Class (TE)",
        "matched_gene_id": "Gene", "strand_match": "Strand match",
        "confirmed_by_junction_screen": "Junction-confirmed", "max_tpm": "Max TPM",
    }

    thead = "".join(f"<th>{esc(col_label[c])}</th>" for c in cols)
    trows = []
    for r in top:
        tpm = max_tpm.get(r.get("transcript_id", ""), 0.0)
        cells = []
        for c in cols:
            if c == "max_tpm":
                cells.append(f"<td>{tpm:.2f}</td>")
            elif c == "confirmed_by_junction_screen":
                v = r.get(c, "no")
                mark = "&#9733;" if v == "yes" else ""
                cells.append(f"<td>{esc(v)} {mark}</td>")
            else:
                cells.append(f"<td>{esc(r.get(c, '.'))}</td>")
        trows.append("<tr>" + "".join(cells) + "</tr>")

    n_total = len(rows)
    n_confirmed = sum(1 for r in rows if r.get("confirmed_by_junction_screen") == "yes") if has_cross_evidence else None
    shown_note = (
        f"Showing the top {len(top)} of {n_total} candidates (ranked: junction-confirmed "
        "first, then strand-matched, then highest expression). Full table: "
        "<code>results/chimera_assembly/candidates_with_junction_evidence.tsv.gz</code>."
        if has_cross_evidence else
        f"Showing the top {len(top)} of {n_total} candidates (ranked: strand-matched first, "
        "then highest expression). Full table: "
        "<code>results/chimera_assembly/candidates.tsv.gz</code>."
    )

    confirmed_line = (
        f"<li><strong>{n_confirmed} of {n_total}</strong> candidates are also found by the "
        "chimera-junction (read-level) screen -- independent evidence from two different "
        "methods, and your highest-confidence starting point. Marked &#9733; below.</li>"
        if has_cross_evidence else
        "<li>chimera.junction is currently disabled, so no cross-confirmation is available "
        "-- enabling it adds an independent evidence source for these same candidates.</li>"
    )

    html = f"""
<p>This screen finds gene-TE chimeras from StringTie assembly structure --
it catches TE-initiated/exonized/terminated transcripts spliced via an
ordinary, canonical intron, which the chimera-junction screen structurally
cannot see (it only flags reads STAR can't explain as one linear
alignment). Neither screen alone is complete; a candidate found by both
carries the strongest evidence.</p>

<p><strong>What to look at first:</strong></p>
<ul>
{confirmed_line}
<li>For <code>te_initiated</code> calls, prefer <code>strand_match: yes</code> -- a mismatch
usually means the "gene" hit is a spurious overlap, not real transcript
connectivity.</li>
<li>Check <code>te_hits_all</code> / <code>gene_hits_all</code> in the full table: more
than one entry means the locus is ambiguous (a multi-copy TE family, or
overlapping gene isoforms) -- the single <code>te_id</code>/<code>matched_gene_id</code>
reported is just the first hit, not necessarily the right one.</li>
<li>Candidates are structural calls from the merged assembly, not per-sample
-- check <code>tpm_matrix.tsv.gz</code> for real expression, and require
support across multiple replicates of a condition before trusting a call
as biological rather than assembly noise.</li>
<li><code>unspliced_te_only</code> entries have no splice evidence at all -- useful
for manual triage, not for a confident call on their own.</li>
</ul>

<table class="table" style="width:100%; font-size: 90%;">
<thead><tr>{thead}</tr></thead>
<tbody>
{"".join(trows) if trows else "<tr><td colspan=" + str(len(cols)) + ">No candidates.</td></tr>"}
</tbody>
</table>
<p style="font-size: 85%; color: #888;">{shown_note}</p>
"""

    highlights_doc = {
        "id": "chimera_assembly_highlights",
        "parent_id": "chimera_transcripts",
        "parent_name": "Gene-TE chimeras: transcript evidence",
        "section_name": "What to look at",
        "description": "How to read the chimera-assembly output, and this run's top candidates.",
        "plot_type": "html",
        "data": html,
    }
    os.makedirs(os.path.dirname(args.out_highlights), exist_ok=True)
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
