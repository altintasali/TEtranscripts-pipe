#!/usr/bin/env bash
# Guard 49: report sections are grouped and ordered as intended
#
# Run on its own:   .tests/guards/49_report_sections_are_grouped_and_ordered_as_intended.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The junction screen's output used to render as THREE separate
# top-level sections (Chimera / Gene-TE chimeras / Chimera evidence)
# because each view was added by a different script. They now share
# one parent_id across four emitters + sample_qc.R, and the two
# screens are named for the evidence they use rather than by
# "chimera"/"junction", which said nothing about the difference.
#
# The ordering is also load-bearing: how to read the evidence comes
# first, whether that evidence is independent second, and each
# screen's own contribution after that. Before this, every screen
# opened with its own ranked table on its own key.
#
# Four things here are load-bearing and silent when broken: a
# parent_id typo in ANY emitter re-splits a group; the heatmaps
# belong to chimera_structure, not to either screen; the strandedness
# move relies on report_section_order's inverted module-level
# semantics ("before: rseqc" renders AFTER it); and subsections sort
# by name unless ordered.
mkdir -p "$T/org/qc"
python3 - "$T/org/qc" <<'PY'
import json, sys
d = sys.argv[1]
docs = [
 # ONE group now. Five separate Chimera [...] groups fragmented the TOC, and
 # MultiQC has no third heading level, so Reads / Assembly / Evidence
 # structure are carried by section NAMES plus report_section_order.
 ("chimera_candidates_table", "chimera", "Chimera", "Candidates"),
 ("chimera_signal_guide", "chimera", "Chimera", "How to weigh this evidence"),
 ("chimera_evidence_composition", "chimera", "Chimera", "Evidence composition"),
 ("sample_evidence_status", "chimera", "Chimera", "Per-sample status"),
 ("junction_highlights", "chimera", "Chimera", "Reads - what this screen sees"),
 ("chimera_junction_qc", "chimera", "Chimera", "Reads - composition by direction"),
 ("chimera_canonical_rate", "chimera", "Chimera", "Reads - splice-motif rate by direction"),
 ("chimera_te_gene_chimeras", "chimera", "Chimera", "Reads - gene-TE subset"),
 ("chimera_reads_sample_qc_pca", "chimera", "Chimera", "Reads - PCA"),
 ("chimera_assembly_highlights", "chimera", "Chimera", "Assembly - what this screen sees"),
 ("chimera_assembly_classes", "chimera", "Chimera", "Assembly - composition by class"),
 ("chimera_assembly_strand_rate", "chimera", "Chimera", "Assembly - strand-match rate by class"),
 ("chimera_assembly_sample_qc_pca", "chimera", "Chimera", "Assembly - PCA"),
 ("chimera_evidence_correlation", "chimera", "Chimera", "Evidence structure - correlation"),
 ("chimera_evidence_candidates", "chimera", "Chimera", "Evidence structure - leaders by dimension"),
 ("strandedness_check", "strandedness_check", "Strandedness check", "Declared vs. inferred"),
 ("evidence_overview", "evidence_overview", "TE analysis", "What this run measured"),
]
for did, pid, pname, sec in docs:
    json.dump({"id": did, "parent_id": pid, "parent_name": pname,
               "section_name": sec, "plot_type": "html",
               "data": f"<p>{sec}</p>"}, open(f"{d}/{did}_mqc.json", "w"))
PY
printf '\n\nThis is PairEnd Data\nFraction of reads failed to determine: 0.05\nFraction of reads explained by "1++,1--,2+-,2-+": 0.03\nFraction of reads explained by "1+-,1-+,2++,2--": 0.92\n' \
  > "$T/org/qc/s_infer_experiment.txt"
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/org/out" -n r "$T/org/qc" > "$T/org/log" 2>&1; then
  echo "ERROR: multiqc failed on the reorganised fixture"
  tail -30 "$T/org/log"; FAIL=1
else
  python3 - "$T/org/out/r.html" <<'PY' || FAIL=1
import re, sys
h = open(sys.argv[1]).read()
tops, seen = re.findall(r"<h2[^>]*>\s*(?:<[^>]+>\s*)*([^<]+?)\s*<", h), []
for t in tops:
    t = t.strip()
    if t and t not in seen:
        seen.append(t)
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
check(seen.count("Chimera") == 1,
      f"chimera must be exactly ONE top-level section; sections were {seen}")
# the five groups this replaced, plus every earlier naming generation
for gone in ("Chimera [Candidates]", "Chimera [Evidence structure]",
             "Chimera [Reads]", "Chimera [Assembly]", "Chimera [Per-sample]",
             "Gene-TE chimeras", "Gene-TE chimeras: candidates",
             "Gene-TE chimeras: read evidence",
             "Gene-TE chimeras: transcript evidence",
             "Gene-TE chimeras: reading the evidence",
             "Evidence: chimeric reads", "Evidence: assembled transcripts",
             "How much to trust the ranking", "How independent is this evidence?",
             "Per-sample evidence status", "Chimera evidence", "Start here"):
    check(gone not in seen, f"stale top-level section {gone!r} still present")

# Section order is load-bearing and MultiQC sorts by NAME unless every entry
# is listed in report_section_order -- so this pins the whole sequence, not
# just the first item.
order = []
for m in re.finditer(r'href="#([a-z0-9_\-]+)"', h):
    if m.group(1) not in order:
        order.append(m.group(1))
expected = ["chimera_candidates_table", "chimera_signal_guide",
            "chimera_evidence_composition", "sample_evidence_status",
            "junction_highlights", "chimera_junction_qc",
            "chimera_canonical_rate", "chimera_te_gene_chimeras",
            "chimera_reads_sample_qc_pca", "chimera_assembly_highlights",
            "chimera_assembly_classes", "chimera_assembly_strand_rate",
            "chimera_assembly_sample_qc_pca", "chimera_evidence_correlation",
            "chimera_evidence_candidates"]
present = [s for s in expected if s in order]
check(present == sorted(present, key=order.index),
      f"chimera sections render out of order: {[s for s in order if s in expected]}")
# the candidate list must lead: it is the entry point the rest explains
if present:
    check(present[0] == "chimera_candidates_table",
          f"Candidates must render first; got {present[0]}")

# strandedness sits with RSeQC (inverted before/after semantics)
if "RSeQC" in seen and "Strandedness check" in seen:
    check(seen.index("Strandedness check") == seen.index("RSeQC") + 1,
          f"strandedness must render right after RSeQC; got {seen}")
else:
    check(False, f"expected both RSeQC and Strandedness check; got {seen}")
# reading guide first inside each group, beating the alphabetical default
order = []
for m in re.finditer(r'href="#([a-z0-9_\-]+)"', h):
    if m.group(1) not in order:
        order.append(m.group(1))
sys.exit(0 if ok else 1)
PY
fi

echo "[guards] done, FAIL=$FAIL"
exit $FAIL

exit $FAIL
