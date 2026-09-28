#!/usr/bin/env python3
"""The report's "TE analysis" section: what this run actually measured, and
which parts of it are independent evidence.

Written because the pipeline exposes several layers that are easy to mistake
for corroborating results when they are not.  TEcount and TElocal are two
views of the *same* counts; star.two_pass is an upstream alignment knob, not
a result; and of the chimera outputs only the junction screen and the
assembly screen are genuinely independent of each other.  A reader who takes
"four sections all mention TEs" as four confirmations draws the wrong
conclusion, so the report says plainly which is which.

Runs as a snakemake `script:` (like config_used_mqc.py) so it can read the
resolved switches straight from params rather than re-deriving config.
Emits a single MultiQC custom-content HTML document; MultiQC's custom_content
module requires the payload under a top-level "data" key.
"""
import json
import os


def _row(cells, muted=False):
    style = ' style="color:#888;"' if muted else ""
    return "<tr" + style + ">" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"


def main():
    p = snakemake.params
    telocal = bool(p._telocal_enabled)
    chimeric_reads = bool(p._chimera_chimeric_reads_enabled)
    assembly = bool(p._chimera_assembly_enabled)
    sj = bool(p._chimera_splice_junctions_enabled)
    two_pass = str(p._two_pass)
    n_samples = int(p._sample_count)

    on = "&#10003;"
    off = "&#8212;"

    # --- what ran, and what kind of thing it is ------------------------
    layers = [
        (
            "TEcount", on,
            "Expression, per TE <strong>subfamily</strong> + per gene",
            "Baseline quantification. Every copy of a subfamily pooled.",
            False,
        ),
        (
            "TElocal", on if telocal else off,
            "Expression, per TE <strong>insertion</strong>",
            "Same reads, same EM assignment, finer resolution -- "
            "<strong>not</strong> independent confirmation of TEcount.",
            not telocal,
        ),
        (
            "Chimera (chimeric reads)", on if chimeric_reads else off,
            "<strong>Evidence</strong>: reads STAR cannot align linearly",
            "Annotation-blind. Finds breakpoints; blind to chimeras spliced "
            "through an ordinary intron.",
            not chimeric_reads,
        ),
        (
            "Chimera (assembly)", on if assembly else off,
            "<strong>Evidence</strong>: transcript structure from StringTie",
            "Annotation-guided. Finds canonically spliced chimeras; blind to "
            "structures no assembler would build.",
            not assembly,
        ),
        (
            "Chimera (SJ)", on if sj else off,
            "<strong>Evidence</strong>: STAR's own splice junctions (SJ.out.tab)",
            "Same blind spot as assembly (ordinary-intron TE splices "
            "invisible to the junction screen), caught at the read-junction "
            "level instead -- no assembly needed. Brand new, unvalidated.",
            not sj,
        ),
        (
            "STAR 2-pass", f"{on} ({two_pass})" if two_pass != "none" else off,
            "Alignment setting, not a result",
            "Improves junction detection for everything above. Changes the "
            "numbers; is not itself an answer.",
            two_pass == "none",
        ),
    ]
    layer_rows = "".join(
        _row([name, state, kind, note], muted=muted)
        for name, state, kind, note, muted in layers
    )

    # --- the one line that actually resolves the confusion --------------
    n_screens = sum([chimeric_reads, assembly, sj])
    if n_screens >= 2:
        cross_ref = (
            " Chimeric-reads+assembly agreement additionally has its own "
            "cross-referenced file, <code>candidates_with_junction_evidence."
            "tsv.gz</code>." if chimeric_reads and assembly else ""
        )
        independence = (
            f"<p><strong>{n_screens} independent chimera screens are "
            "running.</strong> They look for gene-TE chimeras in ways that "
            "fail differently, so a candidate found by more than one is "
            "stronger evidence than any single screen alone -- see "
            "<code>results/chimera/candidates.tsv.gz</code>'s found_by/"
            f"evidence columns.{cross_ref} Everything else is single-method "
            "evidence.</p>"
        )
    elif n_screens == 1:
        running = "chimeric reads" if chimeric_reads else "assembly" if assembly else "SJ"
        others = ", ".join(
            f"<code>chimera.{key}.enabled: true</code>"
            for key, on_ in (("chimeric_reads", chimeric_reads), ("assembly", assembly),
                              ("splice_junctions", sj))
            if not on_
        )
        independence = (
            f"<p><strong>One chimera screen is running</strong> ({running}). "
            f"There is no second, independent method to cross-check its "
            f"calls; {others} adds one.</p>"
        )
    else:
        independence = (
            "<p>No chimera screen is enabled -- this run is quantification "
            "only.</p>"
        )

    # --- reading order ---------------------------------------------------
    steps = ["<li>Check the QC first: FastQC and STAR alignment rates, then "
             "<em>Strandedness check</em> - a wrong strandedness call "
             "silently distorts every count below, so it is worth confirming "
             "before reading any of them.</li>",
             "<li>Read expression: <em>TEcount</em> for which subfamilies "
             "move" + (", then <em>TElocal</em> for which copy is "
                       "responsible" if telocal else "") + ".</li>"]
    # BUG FIXED 2026: this used to point at "Chimera -> What to look at" /
    # "Chimera (assembly) -> What to look at", two sections that no longer
    # exist (the report's per-screen sections are now "Chimeric reads - what
    # this screen sees" / "Assembly - what this screen sees", describing each
    # screen's blind spots, not a candidate list -- see guard 38) and called
    # the result "ranked", which the pipeline never does (guard 50). The
    # actual unified, cross-screen, sortable-not-ranked candidate table is
    # the "Candidates" section (chimera_candidates_table_mqc.py); one bullet
    # covers both screens since they share that one table.
    if chimeric_reads or assembly:
        steps.append(
            "<li>Open <em>Candidates</em> for the unified gene-TE junction "
            "table (sortable by evidence count, not ranked) instead of the "
            "raw per-screen catalogs.</li>"
        )
    steps.append(
        "<li>The pipeline does not run differential-expression analysis "
        "itself -- take the per-sample TEcount tables into your own "
        "DESeq2/edgeR analysis downstream.</li>"
    )

    html = f"""
<p>This run processed <strong>{n_samples}</strong> sample(s). It answers two
separate questions -- <em>what is expressed</em>, and <em>where genes and TEs
are fused into one transcript</em> -- using the layers below. They are not
all independent, which is the usual source of confusion:</p>

<table class="table" style="width:100%; font-size: 90%;">
<thead><tr><th>Layer</th><th>Ran</th><th>What it is</th><th>Read it as</th></tr></thead>
<tbody>
{layer_rows}
</tbody>
</table>

{independence}

<p><strong>Suggested reading order:</strong></p>
<ol>
{"".join(steps)}
</ol>

<p style="font-size: 85%; color: #888;">Chimera candidates and the evidence
behind them are in the <strong>Chimera</strong> section below.</p>
"""

    doc = {
        "id": "evidence_overview",
        "parent_id": "evidence_overview",
        "parent_name": "TE analysis",
        "section_name": "What this run measured",
        "description": (
            "The evidence layers in this report, and which of them are "
            "independent of each other."
        ),
        "plot_type": "html",
        "data": html,
    }

    out = snakemake.output[0]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")


# Guarded so the module can be imported (by the unit tests) without
# running. Snakemake's script: directive executes the file with
# __name__ == "__main__", so this still runs under the workflow --
# benchmark_summary.py has been doing exactly this all along.
if __name__ == "__main__":
    main()
