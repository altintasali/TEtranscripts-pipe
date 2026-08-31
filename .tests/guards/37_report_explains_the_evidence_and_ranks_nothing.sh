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
check(g.get("parent_id") == "chimera_evidence_guide",
      f"parent_id must be chimera_evidence_guide, got {g.get('parent_id')!r}")
b = g["data"]
# the stance itself must be stated, not implied
check("does not rank chimera candidates" in b,
      "guide must say plainly that the pipeline does not rank")
check("validate candidates manually" in b,
      "guide must tell the reader validation is theirs")
# every measured finding that justifies the stance must be present
for probe in ("chance rate", "6.7% vs 10.2%", "19,503",
              "Splice motif", "Read depth", "TE locus expressed"):
    check(probe in b, f"guide must render {probe!r}")
# n_evidence must be described as a count, with its known bias
check("not a score" in b, "n_evidence must be labelled a count, not a score")
check("assembly" in b and "two of the four flags" in b,
      "guide must disclose the n_evidence bias toward assembly-found pairs")
# ANTI-REGRESSION: the guide explains signals; it never lists pairs.
for leaked in ("L1PA2_dup1", "AluY_dup9", "MIR_dup2", "Gapdh", "Tp53"):
    check(leaked not in b,
          f"guide must not list candidates; found {leaked!r}")
c = json.load(open(f"{d}/chimera_evidence_composition_mqc.json"))
check(c.get("plot_type") == "bargraph", "composition must be a bargraph")
check(c.get("parent_id") == "chimera_evidence_guide",
      "composition must share the guide group")
counts = c["data"]["Gene-TE pairs"]
check(counts.get("Splice motif") == 2, f"splice-motif count wrong: {counts}")
check(counts.get("Called by both screens") == 1, f"both-screens count wrong: {counts}")
check(counts.get("No evidence flag") == 1, f"no-flag count wrong: {counts}")
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
# An empty table is a normal outcome and must not crash.
printf 'gene_id\tte_id\tevidence\tn_evidence\tfound_by\n' | gzip -c > "$T/cand/empty.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
      --evidence "$T/cand/empty.tsv.gz" \
      --out-guide "$T/cand/empty_mqc.json" \
      --out-composition "$T/cand/empty_comp_mqc.json" > "$T/cand/empty.log" 2>&1; then
  echo "ERROR: chimera_evidence_guide_mqc.py died on an empty table"
  cat "$T/cand/empty.log"; FAIL=1
fi

exit $FAIL
