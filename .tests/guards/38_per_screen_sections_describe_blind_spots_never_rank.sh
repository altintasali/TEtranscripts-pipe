#!/usr/bin/env bash
# Guard 38: per-screen sections describe blind spots, never rank
#
# Run on its own:   .tests/guards/38_per_screen_sections_describe_blind_spots_never_rank.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_evidence

# --- The counterpart to the guard above: each screen keeps only what
# is specific to it. A candidate table reappearing here is the exact
# regression that produced three competing rankings.
mkdir -p "$T/hl/qc"
cp "$T/ev/j.tsv.gz" "$T/hl/te-gene-chimeras.tsv.gz"
if ! python3 workflow/scripts/junction_highlights_mqc.py \
      --te-events "$T/hl/te-gene-chimeras.tsv.gz" \
      --out "$T/hl/qc/junction_highlights_mqc.json" > "$T/hl/log" 2>&1; then
  echo "ERROR: junction_highlights_mqc.py failed"; cat "$T/hl/log"; FAIL=1
elif ! python3 -c "
import json,sys
d=json.load(open('$T/hl/qc/junction_highlights_mqc.json'))
ok=True
if 'data' not in d:
    print('ERROR: custom_content needs a top-level data key'); ok=False
# every read-evidence view shares this group; a mismatch splits the
# section in two. Checked separately so a failure here does not
# print as a content failure (it did, once).
if d.get('parent_id')!='chimera_reads':
    print(f\"ERROR: parent_id must be chimera_reads, got {d.get('parent_id')!r}\"); ok=False
b=d.get('data','')
if 'What it cannot see' not in b:
    print('ERROR: screen section must state its blind spot'); ok=False
if '<table' in b:
    print('ERROR: per-screen section must NOT render a candidate table'); ok=False
if 'MIR_dup2' in b:
    print('ERROR: per-screen section must not list individual candidates'); ok=False
sys.exit(0 if ok else 1)
"; then
  echo "ERROR: junction screen-notes doc is wrong (see above)"; FAIL=1
elif ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/hl/out" -n r "$T/hl/qc" > "$T/hl/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the junction screen notes"
  tail -30 "$T/hl/render.log"; FAIL=1
elif ! grep -q "What this screen sees" "$T/hl/out/r.html"; then
  echo "ERROR: junction screen-notes section not rendered"; FAIL=1
fi
# An empty gene-TE table is a normal outcome and must not crash.
printf 'event_id\tgene_id\tte_id\tcanonical\tn_samples\ttotal_reads\n' \
  | gzip -c > "$T/hl/empty.tsv.gz"
if ! python3 workflow/scripts/junction_highlights_mqc.py \
      --te-events "$T/hl/empty.tsv.gz" --out "$T/hl/empty_mqc.json" \
      > "$T/hl/empty.log" 2>&1; then
  echo "ERROR: junction_highlights_mqc.py died on an empty table"
  cat "$T/hl/empty.log"; FAIL=1
fi

exit $FAIL
