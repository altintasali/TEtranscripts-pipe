#!/usr/bin/env bash
# Guard 36: unified chimera table reports evidence, scores nothing
#
# Run on its own:   .tests/guards/36_unified_chimera_table_reports_evidence_scores_nothing.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_evidence

# --- The unified table is the one place the three screens can be
# compared. It used to carry a four-tier confidence ladder; that was
# removed for lacking any validated weighting, so what is pinned here
# is that it reports evidence and does NOT score: no tier column, read
# depth creates no flag at all (the metric most inflated by artifacts),
# TE expression creates a corroboration flag but never a screen-bound
# one, and the sort is a deterministic count rather than a verdict.
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
# evidence/n_evidence were themselves replaced: they mixed screen-bound
# flags (bounded by n_screens) with cross-cutting ones (not bounded by it),
# which read as a contradiction whenever a single-screen pair out-counted a
# multi-screen one. screen_evidence/corroboration split them.
check("evidence" not in h and "n_evidence" not in h,
      "evidence/n_evidence must not be reintroduced")
check(all(c in h for c in
          ("screen_evidence", "n_screen_evidence", "corroboration", "n_corroboration")),
      "screen_evidence/n_screen_evidence/corroboration/n_corroboration columns missing")
g = d.get(("Gapdh", "L1PA2_dup1"), {})
check(g.get("found_by") == "cr+assembly", "Gapdh pair not marked found_by cr+assembly")
check(set(g.get("screen_evidence", "").split(",")) ==
      {"cr_canonical", "assembly_strand_match"},
      f"Gapdh pair screen_evidence set wrong: {g.get('screen_evidence')!r}")
check(g.get("n_screen_evidence") == "2",
      f"n_screen_evidence must count screen-bound flags, got {g.get('n_screen_evidence')}")
check(set(g.get("corroboration", "").split(",")) ==
      {"multi_sample", "telocal_expressed"},
      f"Gapdh pair corroboration set wrong: {g.get('corroboration')!r}")
check(g.get("n_corroboration") == "2",
      f"n_corroboration must count cross-cutting flags, got {g.get('n_corroboration')}")
check(g.get("n_screens") == "2",
      "Gapdh pair found by cr+assembly -- n_screens must count that "
      "WITHOUT it also inflating either evidence count (the old "
      "both_screens double-count)")
check(g.get("cr_events") == "2", "two junction events must collapse into one pair row")
check(g.get("cr_reads") == "60", "junction reads must sum across events")
check(g.get("cr_canonical") == "yes", "canonical on any event must set the pair canonical")
# TE expression IS corroboration at this stage (not screen evidence: TElocal
# is a fourth data source, not one of the three detection screens). One
# small run suggested it discriminates nothing, but that is not enough to
# demote it -- see the evidence guide. Pinned so the flag is neither
# dropped nor silently renamed.
a = d.get(("Actb", "AluY_dup9"), {})
check(a.get("telocal_active") == "yes", "fixture should have an expressed locus here")
check(a.get("corroboration") == "multi_sample,telocal_expressed",
      f"an expressed locus must add the telocal flag; got {a.get('corroboration')!r}")
check(a.get("screen_evidence") == ".",
      f"Actb pair has no canonical motif and no assembly hit; got {a.get('screen_evidence')!r}")
check(d.get(("Myc", "L1MdA_dup4"), {}).get("screen_evidence") == "cr_canonical",
      "canonical-only pair must carry exactly that screen-evidence flag")
check(d.get(("Myc", "L1MdA_dup4"), {}).get("corroboration") == ".",
      "canonical-only pair must carry no corroboration")
# depth is the metric artifacts inflate most: 999 reads earns nothing
tp = d.get(("Tp53", "MIR_dup2"), {})
check(tp.get("screen_evidence") == ".", "999 reads must produce NO screen-evidence flag")
check(tp.get("n_screen_evidence") == "0", "999 reads must leave n_screen_evidence at 0")
check(tp.get("corroboration") == ".", "999 reads must produce NO corroboration flag")
check(tp.get("n_corroboration") == "0", "999 reads must leave n_corroboration at 0")
s = d.get(("Sox2", "SVA_dup3"), {})
check(s.get("found_by") == "assembly", "assembly-only pair mislabelled")
check(s.get("telocal_active") == ".", "assembly-only pair must not claim a telocal verdict")
check(s.get("screen_evidence") == "assembly_strand_match",
      "assembly-only pair's strand match must still be screen evidence")
# deterministic: (n_screens, n_screen_evidence, n_corroboration) desc, then gene
order = [(int(r[h.index("n_screens")]), int(r[h.index("n_screen_evidence")]),
          int(r[h.index("n_corroboration")]), r[h.index("gene_id")]) for r in rows]
check(order == sorted(order, key=lambda x: (-x[0], -x[1], -x[2], x[3])),
      "rows must sort by n_screens desc, n_screen_evidence desc, "
      f"n_corroboration desc, then gene_id; got {order}")
sys.exit(0 if ok else 1)
PY
fi
# ...and it must still work with the assembly screen off.
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/ev/j.tsv.gz" \
      --out "$T/ev/j_only.tsv.gz" > "$T/ev/log2" 2>&1; then
  echo "ERROR: chimera_evidence.py failed without --assembly"; cat "$T/ev/log2"; FAIL=1
elif gzip -dc "$T/ev/j_only.tsv.gz" \
     | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="found_by") c=i; next} {print $c}' \
     | grep -qv '^cr$'; then
  echo "ERROR: without --assembly every pair must be found_by cr"; FAIL=1
fi

exit $FAIL
