#!/usr/bin/env bash
# Guard 36: unified chimera table reports evidence, scores nothing
#
# Run on its own:   .tests/guards/36_unified_chimera_table_reports_evidence_scores_nothing.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_evidence

# --- The unified table is the one place the two screens can be
# compared. It used to carry a four-tier confidence ladder; that was
# removed for lacking any validated weighting, so what is pinned here
# is that it reports evidence and does NOT score: no tier column, the
# two non-evidence signals (read depth, TE expression) create no
# flag, and the sort is a deterministic count rather than a verdict.
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/ev/j.tsv.gz" \
      --assembly "$T/ev/a.tsv.gz" --out "$T/ev/out.tsv.gz" > "$T/ev/log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev/log"; FAIL=1
else
  python3 - "$T/ev/out.tsv.gz" <<'PY' || FAIL=1
import gzip, sys
rows = [l.rstrip("\n").split("\t") for l in gzip.open(sys.argv[1], "rt")]
h, rows = rows[0], rows[1:]
d = {(r[h.index("gene_id")], r[h.index("te_id")]): dict(zip(h, r)) for r in rows}
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
check(len(rows) == 5, f"expected 5 gene-TE pairs, got {len(rows)}")
# the scoring column must be GONE, not merely unused
check("confidence_tier" not in h, "confidence_tier must not be reintroduced")
check("evidence" in h and "n_evidence" in h, "evidence/n_evidence columns missing")
g = d.get(("Gapdh", "L1PA2_dup1"), {})
check(g.get("found_by") == "both", "Gapdh pair not marked found_by both")
check(set(g.get("evidence", "").split(",")) ==
      {"canonical", "multi_sample", "both_screens", "assembly_strand_match",
       "telocal_expressed"},
      f"Gapdh pair evidence set wrong: {g.get('evidence')!r}")
check(g.get("n_evidence") == "5", f"n_evidence must count flags, got {g.get('n_evidence')}")
check(g.get("junction_events") == "2", "two junction events must collapse into one pair row")
check(g.get("junction_reads") == "60", "junction reads must sum across events")
check(g.get("junction_canonical") == "yes", "canonical on any event must set the pair canonical")
# TE expression IS an evidence flag at this stage. One small run suggested it
# discriminates nothing, but that is not enough to demote it -- see the
# evidence guide. Pinned so the flag is neither dropped nor silently renamed.
a = d.get(("Actb", "AluY_dup9"), {})
check(a.get("telocal_active") == "yes", "fixture should have an expressed locus here")
check(a.get("evidence") == "multi_sample,telocal_expressed",
      f"an expressed locus must add the telocal flag; got {a.get('evidence')!r}")
check(d.get(("Myc", "L1MdA_dup4"), {}).get("evidence") == "canonical",
      "canonical-only pair must carry exactly that flag")
# depth is the metric artifacts inflate most: 999 reads earns nothing
tp = d.get(("Tp53", "MIR_dup2"), {})
check(tp.get("evidence") == ".", "999 reads must produce NO evidence flag")
check(tp.get("n_evidence") == "0", "999 reads must leave n_evidence at 0")
s = d.get(("Sox2", "SVA_dup3"), {})
check(s.get("found_by") == "assembly", "assembly-only pair mislabelled")
check(s.get("telocal_active") == ".", "assembly-only pair must not claim a telocal verdict")
# deterministic: n_evidence desc, then gene, then te
order = [(int(r[h.index("n_evidence")]), r[h.index("gene_id")]) for r in rows]
check(order == sorted(order, key=lambda x: (-x[0], x[1])),
      f"rows must sort by n_evidence desc then gene_id; got {order}")
sys.exit(0 if ok else 1)
PY
fi
# ...and it must still work with the assembly screen off.
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/ev/j.tsv.gz" \
      --out "$T/ev/j_only.tsv.gz" > "$T/ev/log2" 2>&1; then
  echo "ERROR: chimera_evidence.py failed without --assembly"; cat "$T/ev/log2"; FAIL=1
elif gzip -dc "$T/ev/j_only.tsv.gz" | tail -n +2 | cut -f6 | grep -qv '^junction$'; then
  echo "ERROR: without --assembly every pair must be found_by junction"; FAIL=1
fi

exit $FAIL
