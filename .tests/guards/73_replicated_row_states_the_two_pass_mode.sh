#!/usr/bin/env bash
# Guard 73: the guide's Replicated row states this run's STAR 2-pass mode
#
# Under star.two_pass: cohort, junctions from every sample's first pass are
# added to every sample's index, so a junction is easier to detect again in
# the other samples: SJ samples and Replicated are not fully independent
# detections. Measured on a real 4-sample run against per_sample 2-pass
# (docs/chimera-evidence.md): small, but only ever upward. The guide's
# Replicated row now says so for cohort runs, says the opposite for
# per_sample / none runs, and keeps a generic caveat when the mode is not
# passed -- and the rule passes the configured mode.
#
# Run on its own:   .tests/guards/73_replicated_row_states_the_two_pass_mode.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

if ! sed -n '/^rule chimera_evidence_guide:/,/^rule /p' workflow/rules/chimera_chimeric_reads.smk \
     | grep -q -- '--star-two-pass {params.star_two_pass}'; then
  echo "ERROR: rule chimera_evidence_guide must pass star.two_pass to the guide"
  FAIL=1
fi

fixture_evidence
for mode in cohort per_sample none generic; do
  extra=()
  [ "$mode" != generic ] && extra=(--star-two-pass "$mode")
  if ! python3 workflow/scripts/chimera_evidence_guide_mqc.py \
        --evidence "$T/ev/out.tsv.gz" --sj-require-canonical true "${extra[@]}" \
        --out-guide "$T/guide_$mode.json" --out-composition "$T/comp_$mode.json" \
        > "$T/$mode.log" 2>&1; then
    echo "ERROR: chimera_evidence_guide_mqc.py ($mode) failed"; cat "$T/$mode.log"; exit 1
  fi
done

python3 - "$T" <<'PY' || FAIL=1
import json, re, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
def replicated_row(mode):
    h = json.load(open(f"{T}/guide_{mode}.json"))["helptext"]
    rows = re.findall(r"<tr>.*?</tr>", h, re.S)
    row = [r for r in rows if "Replicated" in r]
    check(len(row) == 1, f"{mode}: expected one Replicated row, got {len(row)}")
    return row[0] if row else ""
co = replicated_row("cohort")
check("This run used <code>star.two_pass: cohort</code>" in co
      and "not fully independent" in co,
      "cohort: the Replicated row must say this run's SJ samples / Replicated "
      "are not fully independent detections")
for mode in ("per_sample", "none"):
    r = replicated_row(mode)
    check(f"This run used <code>star.two_pass: {mode}</code>" in r
          and "independently" in r and "not fully independent" not in r,
          f"{mode}: the Replicated row must say junctions were detected independently")
g = replicated_row("generic")
check("not fully independent" in g and "This run used" not in g,
      "no mode given: a generic cohort caveat, not a claim about this run")
sys.exit(0 if ok else 1)
PY

exit $FAIL
