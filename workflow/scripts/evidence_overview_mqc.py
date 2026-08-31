#!/usr/bin/env python3
"""The report's "Start here" section: what this run actually measured, and
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
    junction = bool(p._chimera_junction_enabled)
    assembly = bool(p._chimera_assembly_enabled)
    two_pass = str(p._two_pass)
    n_samples = int(p._sample_count)
    has_condition = bool(p._has_condition)

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
            "Chimera (junction)", on if junction else off,
            "<strong>Evidence</strong>: reads STAR cannot align linearly",
            "Annotation-blind. Finds breakpoints; blind to chimeras spliced "
            "through an ordinary intron.",
            not junction,
        ),
        (
            "Chimera (assembly)", on if assembly else off,
            "<strong>Evidence</strong>: transcript structure from StringTie",
            "Annotation-guided. Finds canonically spliced chimeras; blind to "
            "structures no assembler would build.",
            not assembly,
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
    if junction and assembly:
        independence = (
            "<p><strong>Two independent screens are running.</strong> They "
            "look for gene-TE chimeras in ways that fail differently, so a "
            "candidate found by <em>both</em> is the strongest call this "
            "pipeline makes -- see <code>candidates_with_junction_evidence."
            "tsv.gz</code>. Everything else is single-method evidence.</p>"
        )
    elif junction:
        independence = (
            "<p><strong>One chimera screen is running</strong> (junction). "
            "There is no second, independent method to cross-check its calls; "
            "<code>chimera.assembly.enabled: true</code> adds one.</p>"
        )
    elif assembly:
        independence = (
            "<p><strong>One chimera screen is running</strong> (assembly). "
            "There is no second, independent method to cross-check its calls; "
            "<code>chimera.junction.enabled: true</code> adds one.</p>"
        )
    else:
        independence = (
            "<p>No chimera screen is enabled -- this run is quantification "
            "only.</p>"
        )

    # --- reading order ---------------------------------------------------
    steps = ["<li>Check the QC first: FastQC and STAR alignment rates, then "
             "<em>Strandedness check</em> \u2014 a wrong strandedness call "
             "silently distorts every count below, so it is worth confirming "
             "before reading any of them.</li>",
             "<li>Read expression: <em>TEcount</em> for which subfamilies "
             "move" + (", then <em>TElocal</em> for which copy is "
                       "responsible" if telocal else "") + ".</li>"]
    if junction:
        steps.append(
            "<li>Open <em>Chimera &rarr; What to look at</em> for the ranked "
            "gene-TE junctions, not the raw catalog.</li>"
        )
    if assembly:
        steps.append(
            "<li>Open <em>Chimera (assembly) &rarr; What to look at</em>; "
            "prefer candidates marked as junction-confirmed.</li>"
        )
    if has_condition:
        steps.append(
            "<li>Differential results are in "
            "<code>results/tetranscripts/</code> -- the report does not "
            "render them.</li>"
        )
    else:
        steps.append(
            "<li>No <code>condition</code> column in the sample sheet, so no "
            "differential comparison was run.</li>"
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

<p style="font-size: 85%; color: #888;">This pipeline reports chimera evidence
but does not rank or score candidates &mdash; see <strong>Gene-TE chimeras:
reading the evidence</strong> below for what each signal is worth, and expect
to validate calls manually.</p>
"""

    doc = {
        "id": "evidence_overview",
        "parent_id": "evidence_overview",
        "parent_name": "Start here",
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
