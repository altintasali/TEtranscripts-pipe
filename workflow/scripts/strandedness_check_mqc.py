#!/usr/bin/env python3
"""Strandedness check: what you declared vs. what RSeQC infers, per sample.

A wrong strandedness is the most damaging silent error in the pipeline. It
never fails a job -- it just makes TEcount count the wrong strand, which
roughly halves gene counts and inverts the antisense signal, and every
downstream number inherits that. It is also easy to get wrong: kit
documentation says "stranded" without saying which direction, and
"forward"/"reverse" mean opposite things in different vocabularies.

RSeQC already infers strandedness for samples set to "auto". This section
compares that inference against a value the sample sheet declared, for every
sample RSeQC ran on (all of them under strandedness.check_provided, the
default), and flags disagreements the way nf-core/rnaseq's strandedness
check does.

The inference is REPORTED, never applied: an explicit sample sheet value
still decides what TEcount is run with. This section tells you to go look;
it does not quietly change your analysis under you.

Emits a MultiQC custom-content table. Note the RSeQC fractions are reported
as the pipeline's own vocabulary (forward/reverse), not RSeQC's raw
"1++,1--,2+-,2-+" pattern strings, so the table can be read without knowing
the encoding.
"""
import json
import os
import re

# Sample sheet vocabulary is normalised to TEtranscripts' own before it
# reaches here ("unstranded" -> "no"); display the pipeline's vocabulary so
# this table matches the --stranded value actually passed to TEcount.
DISPLAY = {"no": "unstranded", "forward": "forward", "reverse": "reverse"}


def parse_fractions(path):
    """forward/reverse/failed fractions from an infer_experiment.py report.

    RSeQC labels its two strand patterns differently for single- and
    paired-end data ("++,--" vs "1++,1--,2+-,2-+"), so key off the first
    token's suffix rather than the whole pattern -- the same rule
    determine_strandedness.py uses, kept deliberately identical so the two
    can never disagree about what the report says.
    """
    forward = reverse = failed = None
    with open(path) as fh:
        for line in fh:
            m = re.search(r"failed to determine:\s*([\d.]+)", line)
            if m:
                failed = float(m.group(1))
                continue
            m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
            if not m:
                continue
            first_token = m.group(1).split(",")[0]  # "++" / "1++" / "+-" / "1+-"
            if first_token.endswith("++"):
                forward = float(m.group(2))
            elif first_token.endswith("+-"):
                reverse = float(m.group(2))
    return forward, reverse, failed


def main():
    p = snakemake.params
    samples = list(p.samples)
    declared = dict(p.declared)  # sample -> "auto"/"no"/"forward"/"reverse"
    reports = dict(zip(samples, snakemake.input.reports))
    calls = dict(zip(samples, snakemake.input.calls))

    data = {}
    n_mismatch = 0
    n_undetermined = 0
    for sample in samples:
        forward, reverse, failed = parse_fractions(reports[sample])
        with open(calls[sample]) as fh:
            inferred = fh.read().strip() or "no"

        want = declared.get(sample, "auto")
        if want == "auto":
            status, used = "auto-detected", inferred
        elif want == inferred:
            status, used = "OK", want
        else:
            status, used = "MISMATCH", want
            n_mismatch += 1

        # An unstranded call can mean "genuinely unstranded" or "RSeQC could
        # not tell" -- the distinction matters when reading a MISMATCH, so
        # surface it rather than leaving both as "unstranded".
        if inferred == "no" and forward is not None and reverse is not None:
            if max(forward, reverse) < 0.6 and abs(forward - reverse) < 0.2:
                n_undetermined += 1

        data[sample] = {
            "declared": DISPLAY.get(want, want),
            "inferred": DISPLAY.get(inferred, inferred),
            "used": DISPLAY.get(used, used),
            "status": status,
            "frac_forward": forward if forward is not None else 0.0,
            "frac_reverse": reverse if reverse is not None else 0.0,
            "frac_failed": failed if failed is not None else 0.0,
        }

    if n_mismatch:
        verdict = (
            f"<strong style='color:#a94442;'>{n_mismatch} of {len(samples)} "
            "sample(s) disagree with the sample sheet.</strong> The declared "
            "value is still what was used — nothing was overridden — "
            "so if the inference is right, these counts are wrong and the run "
            "needs repeating with a corrected samples.csv. Check the kit "
            "documentation before assuming RSeQC is wrong: a low "
            "<em>failed</em> fraction with one strand clearly dominant is "
            "strong evidence."
        )
    else:
        verdict = (
            f"All {len(samples)} sample(s) agree with the sample sheet (or "
            "were auto-detected)."
        )
    if n_undetermined:
        verdict += (
            f" {n_undetermined} sample(s) had no clear strand signal at all "
            "(both fractions low and close) — those were called "
            "unstranded because nothing else could be justified, which is "
            "different from being confidently unstranded."
        )

    doc = {
        "id": "strandedness_check",
        "parent_id": "strandedness_check",
        "parent_name": "Strandedness check",
        "section_name": "Declared vs. inferred",
        "description": (
            "What each sample's strandedness was set to in samples.csv, what "
            "RSeQC's infer_experiment.py says the reads actually look like, "
            "and which value was used for quantification. "
            f"<br><br>{verdict}"
            "<br><br><em>Why this matters:</em> strandedness never fails a "
            "job, it just changes the answer. Counting the wrong strand "
            "roughly halves gene counts and inverts the antisense signal, and "
            "everything downstream — TE quantification, the chimera "
            "screens' strand-match tests, differential results — "
            "inherits that silently."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "strandedness_check_table",
            "title": "Strandedness check",
            "col1_header": "Sample",
        },
        "headers": {
            "declared": {
                "title": "Declared",
                "description": "samples.csv 'strandedness' column "
                               "('auto' when left blank)",
            },
            "inferred": {
                "title": "Inferred",
                "description": "RSeQC infer_experiment.py call",
            },
            "used": {
                "title": "Used",
                "description": "The --stranded value actually passed to "
                               "TEcount/TElocal",
            },
            "status": {
                "title": "Status",
                "description": "MISMATCH = the sample sheet and the data "
                               "disagree; the sample sheet was used",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "OK"}],
                    "warn": [{"s_eq": "auto-detected"}],
                    "fail": [{"s_eq": "MISMATCH"}],
                },
                "cond_formatting_colours": [
                    {"pass": "#5cb85c"},
                    {"warn": "#f0ad4e"},
                    {"fail": "#d9534f"},
                ],
            },
            "frac_forward": {
                "title": "Forward frac",
                "description": "Fraction of reads consistent with a forward "
                               "-stranded library",
                "format": "{:.3f}", "min": 0, "max": 1,
            },
            "frac_reverse": {
                "title": "Reverse frac",
                "description": "Fraction of reads consistent with a reverse "
                               "-stranded library",
                "format": "{:.3f}", "min": 0, "max": 1,
            },
            "frac_failed": {
                "title": "Failed frac",
                "description": "Fraction RSeQC could not assign to either "
                               "(high values undermine the call)",
                "format": "{:.3f}", "min": 0, "max": 1,
            },
        },
        "data": data,
    }

    out = snakemake.output[0]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")

    print(
        f"strandedness check: {len(samples)} sample(s), "
        f"{n_mismatch} mismatch(es), {n_undetermined} without a clear signal"
    )


# Guarded so the module can be imported (by the unit tests) without
# running. Snakemake's script: directive executes the file with
# __name__ == "__main__", so this still runs under the workflow --
# benchmark_summary.py has been doing exactly this all along.
if __name__ == "__main__":
    main()
