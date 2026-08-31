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
 ("chimera_evidence_guide", "chimera_evidence_guide", "Chimera [Candidates]", "How to weigh this evidence"),
 ("chimera_evidence_composition", "chimera_evidence_guide", "Chimera [Candidates]", "Evidence composition"),
 ("chimera_evidence_correlation", "chimera_structure", "Chimera [Evidence structure]", "Evidence correlation"),
 ("chimera_evidence_candidates", "chimera_structure", "Chimera [Evidence structure]", "Candidates by evidence type"),
 ("junction_highlights", "chimera_reads", "Chimera [Reads]", "What this screen sees"),
 ("chimera_junction_qc", "chimera_reads", "Chimera [Reads]", "Events by direction"),
 ("chimera_canonical_rate", "chimera_reads", "Chimera [Reads]", "Splice-motif rate by direction"),
 ("chimera_te_gene_chimeras", "chimera_reads", "Chimera [Reads]", "Gene-TE subset"),
 ("chimera_reads_sample_qc_pca", "chimera_reads", "Chimera [Reads]", "PCA"),
 ("chimera_assembly_highlights", "chimera_transcripts", "Chimera [Assembly]", "What this screen sees"),
 ("chimera_assembly_classes", "chimera_transcripts", "Chimera [Assembly]", "Candidates by class"),
 ("sample_evidence_status", "evidence_status", "Chimera [Per-sample]", "Per-sample status"),
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
check(seen.count("Chimera [Reads]") == 1,
      f"read-evidence must be exactly ONE section; sections were {seen}")
check("Chimera [Assembly]" in seen,
      "transcript-evidence section missing")
check("Chimera [Candidates]" in seen,
      "evidence-guide section missing")
check("Chimera [Evidence structure]" in seen,
      "evidence-independence section missing")
# how to read the evidence, then whether it is independent, then what
# each screen contributed -- the whole point of the ordering.
if all(s in seen for s in ("Chimera [Candidates]",
                           "Chimera [Evidence structure]",
                           "Chimera [Reads]")):
    check(seen.index("Chimera [Candidates]")
          < seen.index("Chimera [Evidence structure]")
          < seen.index("Chimera [Reads]"),
          f"guide -> independence -> screens order broken; got {seen}")
# old headings must be gone as TOP-LEVEL sections. "candidates" and
# "trust the ranking" are here because they shipped briefly and both
# implied a verdict the pipeline no longer issues.
for gone in ("Chimera", "Gene-TE chimeras", "Chimera evidence",
             "Chimera (assembly)", "Gene-TE chimeras: read evidence",
             "Gene-TE chimeras: transcript evidence",
             "Gene-TE chimeras: candidates",
             "How much to trust the ranking"):
    check(gone not in seen, f"stale top-level section {gone!r} still present")
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
for first, later in (("junction_highlights", "chimera_junction_qc"),
                     ("chimera_evidence_guide", "chimera_evidence_composition"),
                     ("chimera_assembly_highlights", "chimera_assembly_classes")):
    if first in order and later in order:
        check(order.index(first) < order.index(later),
              f"{first} must render before {later}")
    else:
        check(False, f"missing {first}/{later} in the report")
sys.exit(0 if ok else 1)
PY
fi

echo "[guards] done, FAIL=$FAIL"
exit $FAIL

exit $FAIL
