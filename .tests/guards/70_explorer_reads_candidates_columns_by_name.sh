#!/usr/bin/env bash
# Guard 70: the explorer rule reads candidates.tsv.gz columns by NAME
#
# chimera_candidates_explorer's shell step pulls gene_id, te_id,
# telocal_locus and assembly_transcript_ids out of candidates.tsv.gz to join
# genome loci and cohort read totals. It used fixed `cut -fN` positions,
# which went stale silently as chimera_evidence.py's OUT_COLUMNS grew: on a
# real run the "TElocal keys" were chimera types (te_terminated) and the
# "assembly keys" were counts (0), so the explorer's TElocal reads and
# Assembly reads were blank for every pair. This runs the rule's own
# _evidence_column() (extracted from the .smk, not a copy) on a fixture whose
# columns are deliberately out of the real order.
#
# Run on its own:   .tests/guards/70_explorer_reads_candidates_columns_by_name.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

SMK=workflow/rules/chimera_chimeric_reads.smk
if grep -n 'evidence} | tail -n +2 | cut -f' "$SMK"; then
  echo "ERROR: $SMK reads candidates.tsv.gz by column position again (cut -fN)"
  FAIL=1
fi

{ printf 'assembly_transcript_ids\tte_id\ttelocal_locus\tgene_id\n'
  printf 'MSTRG.1.1,MSTRG.1.2\tL1_dup1\tL1_dup1:L1:L1:LINE\tG1\n'
  printf '.\tB1_dup2\t.\tG2\n'
} | gzip -c > "$T/cand.tsv.gz"

python3 - "$T" "$SMK" <<'PY' || FAIL=1
import re, subprocess, sys, types
T, smk = sys.argv[1], sys.argv[2]
src = open(smk).read()
m = re.search(r"^def _evidence_column\(name\):\n(?:(?:    .*|)\n)+", src, re.M)
if not m:
    print("ERROR: _evidence_column() not found in", smk); sys.exit(1)
ns = {}
exec(m.group(0), ns)
inp = types.SimpleNamespace(evidence=f"{T}/cand.tsv.gz")
ok = True
def run(col):
    cmd = ns["_evidence_column"](col).format(input=inp)
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                          capture_output=True, text=True)
for col, want in (("gene_id", ["G1", "G2"]), ("te_id", ["L1_dup1", "B1_dup2"]),
                  ("telocal_locus", ["L1_dup1:L1:L1:LINE", "."]),
                  ("assembly_transcript_ids", ["MSTRG.1.1,MSTRG.1.2", "."])):
    r = run(col)
    got = r.stdout.split("\n")[:-1]
    if r.returncode != 0 or got != want:
        print(f"ERROR: column {col!r}: want {want}, got {got} (exit {r.returncode}; {r.stderr.strip()})")
        ok = False
r = run("no_such_column")
if r.returncode == 0:
    print("ERROR: a missing column must fail the rule, not print whole rows")
    ok = False
sys.exit(0 if ok else 1)
PY

exit $FAIL
