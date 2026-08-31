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
# states what is known about each and stops. Pinned below: the
# measured caveats actually render, the composition counts are real,
# and -- the anti-regression that matters -- the guide lists NO
# candidates, because any ordering of rows reads as importance.
mkdir -p "$T/cand/qc"
gzip -dc "$T/ev/out.tsv.gz" | gzip -c > "$T/cand/evidence.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
      --evidence "$T/cand/evidence.tsv.gz" \
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
b = g["data"]
# The "we don't rank" statement now lives ONCE, in the Candidates section --
# it used to be repeated across five sections, which read as nagging next to
# sortable content. The guide's job is explaining what each signal is worth.
check("does not rank" not in b,
      "the no-ranking statement belongs in Candidates, not repeated here")
check("Candidates" in b, "guide must point at the table it explains")
# every measured finding that justifies the stance must still be present
for probe in ("chance rate", "6.7% vs 10.2%", "19,503",
              "Splice motif", "Read depth", "TE locus expressed"):
    check(probe in b, f"guide must render {probe!r}")
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
counts = {k: v["Gene-TE pairs"] for k, v in c["data"].items()}
check(counts.get("Splice motif") == 2, f"splice-motif count wrong: {counts}")
check(counts.get("Called by both screens") == 1, f"both-screens count wrong: {counts}")
check(counts.get("No evidence flag") == 1, f"no-flag count wrong: {counts}")
check("do not sum" in b or "overlap" in b,
      "the section must say the bars overlap rather than partition the cohort")
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
# An empty run is a normal outcome. It must not crash the SCRIPT -- and, the
# part that was missing, it must not crash MULTIQC either: a bargraph whose
# every value is zero raises "No datasets to plot", exits non-zero, and
# produces NO REPORT AT ALL. Rendering the empty case is the only way that
# regression is visible; asserting the script survived is not enough.
mkdir -p "$T/cand/empty_qc"
printf 'gene_id\tte_id\tevidence\tn_evidence\tfound_by\n' \
  | gzip -c > "$T/cand/empty.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
      --evidence "$T/cand/empty.tsv.gz" \
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

exit $FAIL