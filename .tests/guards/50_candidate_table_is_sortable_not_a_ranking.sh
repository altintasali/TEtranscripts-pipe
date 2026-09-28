#!/usr/bin/env bash
# Guard 50: candidate table is a sortable view, not a ranking
#
# Run on its own:   .tests/guards/50_candidate_table_is_sortable_not_a_ranking.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_evidence

# --- The report needs an entry point users can act on, and this pipeline has
# twice had a RANKING removed for having no validated weighting. Both have to
# stay true at once, which is what this pins: the table exists and is usable,
# its columns are grouped by screen (CR / SJ / Assembly / Support) instead of
# interleaved, and its order is a single COUNT the reader can re-sort rather
# than a score the pipeline asserts.
mkdir -p "$T/tbl/qc"
gzip -dc "$T/ev/out.tsv.gz" | gzip -c > "$T/tbl/evidence.tsv.gz"
{ printf 'gene_id\tgene_name\n'; printf 'Gapdh\tGAPDH\n'; } | gzip -c > "$T/tbl/sym.tsv.gz"

if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/tbl/evidence.tsv.gz" --gene-names "$T/tbl/sym.tsv.gz" \
      --out "$T/tbl/qc/chimera_candidates_table_mqc.json" \
      > "$T/tbl/log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed"; cat "$T/tbl/log"; FAIL=1
else
  python3 - "$T/tbl/qc/chimera_candidates_table_mqc.json" <<'PY' || FAIL=1
import json, sys
d = json.load(open(sys.argv[1]))
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False

check(d.get("parent_id") == "chimera",
      f"parent_id must be chimera, got {d.get('parent_id')!r}")
# MultiQC's native table, NOT raw html: only the native table gives the
# reader the sort/filter toolbox, which is the entire premise here.
check(d.get("plot_type") == "table",
      f"must be a MultiQC table so the reader can sort; got {d.get('plot_type')!r}")

# every signal is its own column, grouped Pair / Overview / CR / SJ /
# Assembly / Support, so any of them can be sorted on -- Gene and TE
# insertion included (hidden, but present), so a reader can find a
# specific gene.
headers = d.get("headers", {})
for col in ("Gene", "TE insertion", "TE subfamily", "TE class", "Screens",
            "Found by", "CR motif", "CR samples", "CR reads", "SJ motif",
            "SJ samples", "SJ reads", "Assembly strand",
            "Assembly transcripts", "Replicated", "TElocal reads"):
    check(col in headers, f"column {col!r} missing -- the reader cannot sort on it")

# ANTI-REGRESSION: the two counts this table used to show (each nearly a
# duplicate of information now visible as its own column) must not come
# back as table headers -- they still exist in candidates.tsv.gz and the
# explorer (hidden there), just not rendered here a second time.
for gone in ("Screen evidence count", "Corroboration count"):
    check(gone not in headers, f"removed count column leaked back in: {gone!r}")

# Gene / TE insertion are real columns (sortable/filterable) but hidden from
# the initial view -- both already appear in the row key.
check(headers.get("Gene", {}).get("hidden") is True,
      "Gene column must be hidden by default")
check(headers.get("TE insertion", {}).get("hidden") is True,
      "TE insertion column must be hidden by default")

# the table opens on Screens alone now -- sort_rows does not stick by itself
check(d["pconfig"].get("defaultsort") == [{"column": "Screens", "direction": "desc"}],
      f"defaultsort must be Screens desc only; got {d['pconfig'].get('defaultsort')}")

# rows are sorted by (n_screens desc, gene label asc, te_id asc) -- a sort
# that asserts nothing beyond screen count, unlike candidates.tsv.gz's own
# finer (n_screens, n_screen_evidence, n_corroboration) order.
rows = list(d["data"])
# " | ", never " / ": MultiQC cleans table row names like filenames and splits
# on "/", so a "GENE / te_id" key rendered as just the TE id and the gene
# vanished from the report entirely.
check(all("/" not in r for r in rows),
      f"row keys must not contain '/' -- MultiQC strips everything before it; got {rows[:2]}")
check(rows[0].startswith("GAPDH |"),
      f"gene symbols must be resolved and the only 2-screen pair first; got {rows[0]!r}")
# fixture_evidence: Gapdh(2 screens) first, then the four 1-screen pairs
# alphabetically (Actb, Myc, Sox2, Tp53 -- none has a gene_name, so gene_id
# is the label) -- this is the alphabetical tie-break the rewrite added.
expected_order = ["GAPDH | L1PA2_dup1", "Actb | AluY_dup9", "Myc | L1MdA_dup4",
                   "Sox2 | SVA_dup3", "Tp53 | MIR_dup2"]
check(rows == expected_order,
      f"row order must be (Screens desc, gene asc, TE asc); got {rows}, want {expected_order}")

# Replicated agrees with the multi_sample flag chimera_evidence.py itself
# computed, for every row -- both Gapdh and Actb have a junction with
# n_samples > 1 and must read "yes"; Myc/Tp53/Sox2 must read "no".
want_replicated = {
    "GAPDH | L1PA2_dup1": "yes", "Actb | AluY_dup9": "yes",
    "Myc | L1MdA_dup4": "no", "Sox2 | SVA_dup3": "no", "Tp53 | MIR_dup2": "no",
}
for r, want in want_replicated.items():
    got = d["data"][r].get("Replicated")
    check(got == want, f"Replicated wrong for {r!r}: got {got!r}, want {want!r}")

# TElocal blank-vs-0 survives: Sox2 (assembly-only, telocal never touched
# it) must have NO "TElocal reads" key at all (blank cell), while at least
# one other row (telocal_active resolved to yes/no) has an explicit 0.
check("TElocal reads" not in d["data"]["Sox2 | SVA_dup3"],
      "TElocal reads must be blank (unset), not 0, when TElocal never ran for a pair")
check(any("TElocal reads" in v and v["TElocal reads"] == 0 for k, v in d["data"].items() if k != "Sox2 | SVA_dup3"),
      "expected at least one row with TElocal reads == 0 (measured, found nothing)")

# ANTI-REGRESSION: the ordering must be disclosed as a count, and the section
# must not claim the top rows are correct. A four-key lexicographic sort on
# (canonical, replicates, strand, depth) was removed in d927c8f; this is the
# check that stops it coming back wearing a disclosure line.
desc = d.get("description", "")
check("not a score" in desc, "the section must say the order is not a score")
check("sort" in desc.lower(), "the section must tell the reader they can re-sort")
check("validate candidates manually" in desc,
      "the section must say validation is the reader's job")
for banned in ("ranked by", "confidence tier", "highest confidence"):
    check(banned not in desc.lower(),
          f"section must not present itself as a ranking ({banned!r})")
sys.exit(0 if ok else 1)
PY
fi

# --- ...and it must render, with the table toolbox present.
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/tbl/out" -n r "$T/tbl/qc" > "$T/tbl/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the candidates table"
  tail -30 "$T/tbl/render.log"; FAIL=1
elif ! grep -q "PlotType.TABLE" "$T/tbl/render.log"; then
  echo "ERROR: rendered as something other than a MultiQC table"; FAIL=1
elif ! grep -q "Plot Table Data" "$T/tbl/out/r.html"; then
  echo "ERROR: the table toolbox (sort/filter) is absent from the report"; FAIL=1
elif ! python3 - "$T/tbl/out/r_data/multiqc_chimera_candidates_table_plot.txt" <<'PY2'
import sys
# The bug this catches: MultiQC cleans table row names like filenames. A key
# of "GENE / te_id" was split on "/" and only the basename kept, so every gene
# silently disappeared from the rendered report while the emitted JSON looked
# correct. Only a real render shows it -- assert on MultiQC's OUTPUT.
rows = [l.split("\t") for l in open(sys.argv[1]).read().splitlines()]
head, body = rows[0], rows[1:]
ok = True
if "Gene" not in head:
    print("ERROR: rendered table has no Gene column"); ok = False
gi = head.index("Gene") if "Gene" in head else None
names = [r[0] for r in body]
if not any("GAPDH" in n for n in names):
    print(f"ERROR: gene lost from rendered row names: {names}"); ok = False
if gi is not None and not any(r[gi] == "GAPDH" for r in body):
    print("ERROR: gene lost from the rendered Gene column"); ok = False
if "Screen evidence count" in head or "Corroboration count" in head:
    print(f"ERROR: removed count column rendered in output: {head}"); ok = False
sys.exit(0 if ok else 1)
PY2
then
  echo "ERROR: gene names did not survive MultiQC rendering"; FAIL=1
fi

# --- Top-N tie disclosure: when more pairs share the single highest Screens
# value than fit in top_n, the description must say so, since which of the
# tied pairs get shown is then an alphabetical accident.
mkdir -p "$T/tie/qc"
cols="gene_id\tte_id\tte_subfamily\tte_family\tte_class\tfound_by\tn_screens\tscreen_evidence\tn_screen_evidence\tcorroboration\tn_corroboration\tcr_events\tcr_reads\tcr_max_samples\tcr_canonical\tcr_chimera_types\ttelocal_active\ttelocal_count\ttelocal_locus\tassembly_transcripts\tassembly_chimera_types\tassembly_strand_match\tassembly_transcript_ids\tsj_events\tsj_reads\tsj_max_samples\tsj_canonical\tsj_chimera_types"
{
  printf "%b\n" "$cols"
  # four pairs tied at n_screens=3 (the top value), one at n_screens=1 --
  # with --top-n 2, only Aaa and Bbb (alphabetically first of the four
  # tied) are shown, and Ccc/Ddd/Zzz are not.
  printf 'Aaa\tT1\tsf\tfam\tLINE\tcr+assembly+sj\t3\t.\t0\t.\t0\t1\t10\t1\tno\t.\t.\t0\t.\t1\tte_terminated\t.\tMSTRG.a\t1\t10\t1\tno\t.\n'
  printf 'Bbb\tT2\tsf\tfam\tLINE\tcr+assembly+sj\t3\t.\t0\t.\t0\t1\t10\t1\tno\t.\t.\t0\t.\t1\tte_terminated\t.\tMSTRG.b\t1\t10\t1\tno\t.\n'
  printf 'Ccc\tT3\tsf\tfam\tLINE\tcr+assembly+sj\t3\t.\t0\t.\t0\t1\t10\t1\tno\t.\t.\t0\t.\t1\tte_terminated\t.\tMSTRG.c\t1\t10\t1\tno\t.\n'
  printf 'Ddd\tT4\tsf\tfam\tLINE\tcr+assembly+sj\t3\t.\t0\t.\t0\t1\t10\t1\tno\t.\t.\t0\t.\t1\tte_terminated\t.\tMSTRG.d\t1\t10\t1\tno\t.\n'
  printf 'Zzz\tT5\tsf\tfam\tLINE\tcr\t1\t.\t0\t.\t0\t1\t10\t1\tno\t.\t.\t0\t.\t0\t.\t.\t.\t0\t0\t0\tno\t.\n'
} | gzip -c > "$T/tie/candidates.tsv.gz"

if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/tie/candidates.tsv.gz" --top-n 2 \
      --out "$T/tie/qc/chimera_candidates_table_mqc.json" \
      > "$T/tie/log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed on the tie fixture"
  cat "$T/tie/log"; FAIL=1
else
  python3 - "$T/tie/qc/chimera_candidates_table_mqc.json" <<'PY3' || FAIL=1
import json, sys
d = json.load(open(sys.argv[1]))
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
rows = list(d["data"])
check(rows == ["Aaa | T1", "Bbb | T2"],
      f"expected the alphabetically first 2 of the 4 tied pairs; got {rows}")
desc = d.get("description", "")
check("4 pairs share the top Screens value (3)" in desc,
      f"description must disclose the tie count; got {desc!r}")
check("only the alphabetically first 2" in desc,
      f"description must say how many of the tied pairs are actually shown; got {desc!r}")
sys.exit(0 if ok else 1)
PY3
fi

exit $FAIL
