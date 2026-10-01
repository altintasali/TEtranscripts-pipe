#!/usr/bin/env bash
# Guard 37: report explains the evidence and ranks nothing
#
# Run on its own:   .tests/guards/37_report_explains_the_evidence_and_ranks_nothing.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_evidence

# --- The report used to answer "which chimeras are real?" three
# ways, then briefly with one confidence ladder. Both are gone: no
# weighting of these signals has been validated here, so the report
# states what is known about each and stops. Pinned below: the guide is
# now GENERAL (no run-specific numbers baked into every user's report --
# those moved to docs/chimera-evidence.md), the composition counts are
# real, its bars render in the same order as the Candidates table's
# columns, and -- the anti-regression that matters -- the guide lists NO
# candidates, because any ordering of rows reads as importance.
mkdir -p "$T/cand/qc"
gzip -dc "$T/ev/out.tsv.gz" | gzip -c > "$T/cand/evidence.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
      --evidence "$T/cand/evidence.tsv.gz" \
      --sj-require-canonical true \
      --out-guide "$T/cand/qc/chimera_evidence_guide_mqc.json" \
      --out-composition "$T/cand/qc/chimera_evidence_composition_mqc.json" \
      > "$T/cand/log" 2>&1; then
  echo "ERROR: chimera_evidence_guide_mqc.py failed"; cat "$T/cand/log"; FAIL=1
else
  python3 - "$T/cand/qc" <<'PY2' || FAIL=1
import json, sys
d = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
g = json.load(open(f"{d}/chimera_evidence_guide_mqc.json"))
check("data" in g, "custom_content needs a top-level data key")
check(g.get("parent_id") == "chimera",
      f"parent_id must be chimera, got {g.get('parent_id')!r}")
# BUG FIXED 2026 (Task 3, report-bloat trim): the guide's detailed
# signal-by-signal table moved from "data" (always visible) into
# "helptext" (collapsed behind a Help button) -- the whole content is still
# there, just not always-rendered, so these checks now read both fields
# together rather than "data" alone.
b = g["data"] + g.get("helptext", "")
# The "we don't rank" statement now lives ONCE, in the Candidates section --
# it used to be repeated across five sections, which read as nagging next to
# sortable content. The guide's job is explaining what each signal is worth.
check("does not rank" not in b,
      "the no-ranking statement belongs in Candidates, not repeated here")
check("Candidates" in b, "guide must point at the table it explains")
check("docs/chimera-evidence.md" in b,
      "guide must point readers at the project's own measurements")

# BUG FIXED 2026: this guide used to print ONE earlier 4-sample mouse run's
# own numbers into every run's report, whatever the reader's own species or
# cohort -- some of them (e.g. "brand new, has not been run on real data")
# were already false by the time this was reviewed. Those measurements now
# live in docs/chimera-evidence.md with their context; the report itself
# must never contain them again.
for gone in ("chance rate", "6.7%", "10.2%", "19,503", "mouse run",
             "brand new"):
    check(gone not in b, f"guide must NOT contain the run-specific finding {gone!r} -- see docs/chimera-evidence.md")

# every signal label from the Candidates table must be named, in the
# Candidates table's own column order (TE orientation, Screens/Found by, CR
# motif, SJ motif, Assembly strand, Replicated, TElocal reads; Read depth last since
# it spans both CR/SJ and is not evidence at all).
order_labels = ["TE orientation", "Screens / Found by", "CR motif",
                "CR max anchor", "SJ motif",
                "SJ unique fraction / overhang",
                "Assembly strand", "Replicated", "TElocal reads",
                "Read depth"]
positions = []
for label in order_labels:
    i = b.find(label)
    check(i != -1, f"guide must render the signal label {label!r}")
    positions.append(i)
check(positions == sorted(positions),
      f"guide rows must appear in Candidates-table column order "
      f"{order_labels}; got positions {positions}")

# and every signal must name the tool it came from
for probe in ("STAR (chimeric junctions)", "StringTie (assembly)",
              "STAR + StringTie", "TElocal"):
    check(probe in b, f"guide must name the source {probe!r}")
# ANTI-REGRESSION: the guide explains signals; the Candidates table lists
# pairs. Keeping them apart is what stops the guide drifting back into a
# second, differently-ordered candidate list.
for leaked in ("L1PA2_dup1", "AluY_dup9", "MIR_dup2", "Gapdh", "Tp53"):
    check(leaked not in b, f"guide must not list candidates; found {leaked!r}")
c = json.load(open(f"{d}/chimera_evidence_composition_mqc.json"))
check(c.get("plot_type") == "bargraph", "composition must be a bargraph")
check(c.get("parent_id") == "chimera",
      "composition must share the chimera group")
# One bar PER EVIDENCE TYPE, unstacked. The first shape made every type a
# category of a single stacked bar, which drew them as slices of a whole --
# but a pair can carry several flags, so the counts overlap and a partition
# is exactly the wrong picture.
check(c["pconfig"].get("stacking") == "group",
      "composition must not stack: the counts overlap, they are not a partition")
# MultiQC alphabetises bar categories by default (bargraph.py's own
# sort_samples: True) -- without this explicitly off, the categories would
# NOT render in FLAGS' own order below. Measured; see also the rendered
# order check after multiqc runs.
check(c["pconfig"].get("sort_samples") is False,
      "composition must set sort_samples: False, or MultiQC alphabetises the bars")
counts = {k: v["Gene-TE pairs"] for k, v in c["data"].items()}
check(counts.get("CR motif") == 2, f"CR-motif count wrong: {counts}")
# "Called by both screens" is gone: it duplicated found_by/n_screens and was
# double-counted into the combined evidence count that has itself since
# been split into screen_evidence/corroboration (see chimera_evidence.py's
# module docstring) -- screen agreement is no longer one of the
# composition bars.
check("Called by both screens" not in counts,
      f"removed bar leaked back into composition: {counts}")
check(counts.get("No evidence flag") == 1, f"no-flag count wrong: {counts}")
# the COMPOSITION section's own text -- this used to search the guide's
# text (b), and only passed because an unrelated guide sentence happened to
# contain the word "overlap"
comp_text = c.get("description", "") + c.get("helptext", "")
check("do not sum" in comp_text or "overlap" in comp_text,
      "the section must say the bars overlap rather than partition the cohort")
# Bar category order (the JSON dict's own insertion order) must match the
# Candidates table's column order: CR motif, SJ motif, Assembly strand,
# Replicated, TElocal expressed, then No evidence flag last.
check(list(c["data"].keys()) == ["CR motif", "SJ motif", "Assembly strand",
                                  "Replicated", "TElocal expressed",
                                  "No evidence flag"],
      f"composition bar order wrong: {list(c['data'].keys())}")
sys.exit(0 if ok else 1)
PY2
fi
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/cand/out" -n r "$T/cand/qc" > "$T/cand/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the evidence guide"
  tail -30 "$T/cand/render.log"; FAIL=1
elif ! grep -q "How to weigh this evidence" "$T/cand/out/r.html"; then
  echo "ERROR: evidence guide section not rendered"; FAIL=1
elif ! grep -q "Evidence composition" "$T/cand/out/r.html"; then
  echo "ERROR: evidence composition section not rendered"; FAIL=1
fi
# The JSON dict order is necessary but not sufficient -- MultiQC's own
# rendering is what the ticket asked to verify. For a HORIZONTAL bar plot
# (this one), Plotly draws trace-array index 0 at the BOTTOM of the chart
# and later indices upward, so the array MultiQC actually renders with is
# the REVERSE of the natural top-to-bottom reading order; reversing it back
# must reproduce the same order pinned above.
if ! python3 - "$T/cand/out/r_data/multiqc_data.json" <<'PY3'
import json, sys
d = json.load(open(sys.argv[1]))
plot = d["report_plot_data"].get("chimera_evidence_composition_plot")
ok = True
if plot is None:
    print("ERROR: chimera_evidence_composition_plot missing from the rendered report")
    ok = False
else:
    samples = plot["datasets"][0].get("samples", [])
    want = ["CR motif", "SJ motif", "Assembly strand", "Replicated",
            "TElocal expressed", "No evidence flag"]
    got_top_to_bottom = list(reversed(samples))
    if got_top_to_bottom != want:
        print(f"ERROR: rendered composition bar order wrong -- got (top to "
              f"bottom) {got_top_to_bottom}, want {want}")
        ok = False
sys.exit(0 if ok else 1)
PY3
then
  echo "ERROR: composition bar order did not survive MultiQC rendering"; FAIL=1
fi
# An empty run is a normal outcome. It must not crash the SCRIPT -- and, the
# part that was missing, it must not crash MULTIQC either: a bargraph whose
# every value is zero raises "No datasets to plot", exits non-zero, and
# produces NO REPORT AT ALL. Rendering the empty case is the only way that
# regression is visible; asserting the script survived is not enough.
mkdir -p "$T/cand/empty_qc"
printf 'gene_id\tte_id\tscreen_evidence\tn_screen_evidence\tcorroboration\tn_corroboration\tfound_by\n' \
  | gzip -c > "$T/cand/empty.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
      --evidence "$T/cand/empty.tsv.gz" \
      --sj-require-canonical false \
      --out-guide "$T/cand/empty_qc/chimera_evidence_guide_mqc.json" \
      --out-composition "$T/cand/empty_qc/chimera_evidence_composition_mqc.json" \
      > "$T/cand/empty.log" 2>&1; then
  echo "ERROR: chimera_evidence_guide_mqc.py died on an empty table"
  cat "$T/cand/empty.log"; FAIL=1
elif ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/cand/empty.tsv.gz" \
      --out "$T/cand/empty_qc/chimera_candidates_table_mqc.json" \
      >> "$T/cand/empty.log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py died on an empty table"
  cat "$T/cand/empty.log"; FAIL=1
elif ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/cand/empty_out" -n r "$T/cand/empty_qc" \
      > "$T/cand/empty_render.log" 2>&1 || [ ! -f "$T/cand/empty_out/r.html" ]; then
  echo "ERROR: a chimera run with no candidates must still produce a report"
  tail -20 "$T/cand/empty_render.log"; FAIL=1
fi

# --- docs/chimera-evidence.md must exist and actually carry the findings
# that used to be hard-coded into the guide, so they are not silently lost.
if [ ! -f docs/chimera-evidence.md ]; then
  echo "ERROR: docs/chimera-evidence.md is missing"; FAIL=1
else
  for probe in "chance rate" "91%" "6.7%" "10.2%" "19,503" "4-sample mouse run" \
               "d927c8f" "0d04e43" "Not yet measured" "SJ-screen agreement" \
               "Three-screen agreement" "Condition-aware replication"; do
    if ! grep -qF "$probe" docs/chimera-evidence.md; then
      echo "ERROR: docs/chimera-evidence.md missing expected content: $probe"
      FAIL=1
    fi
  done
fi

exit $FAIL
