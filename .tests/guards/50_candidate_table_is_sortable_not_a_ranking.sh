#!/usr/bin/env bash
# Guard 50: candidate table is a sortable view, not a ranking
#
# Run on its own:   .tests/guards/50_candidate_table_is_sortable_not_a_ranking.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The report needs an entry point users can act on, and this pipeline has
# twice had a RANKING removed for having no validated weighting. Both have to
# stay true at once, which is what this pins: the table exists and is usable,
# its columns are grouped by screen (CR / SJ / Assembly / Support) instead of
# interleaved, and its order is a single COUNT the reader can re-sort rather
# than a score the pipeline asserts.
fixture_evidence
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
for col in ("Gene", "TE insertion", "TE subfamily", "TE class",
            "TE position", "TE orientation", "Distance", "Screens",
            "Found by", "Types", "CR motif", "CR samples", "CR reads", "SJ motif",
            "SJ unique fraction", "SJ max overhang",
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

# This fixture's 5 pairs are far below MIN_ROWS, so every one of them is
# shown regardless of which Screens tier they're in -- fill-down exhausts
# the whole file. Order: (n_screens desc, gene asc, TE asc).
rows = list(d["data"])
# " | ", never " / ": MultiQC cleans table row names like filenames and splits
# on "/", so a "GENE / te_id" key rendered as just the TE id and the gene
# vanished from the report entirely.
check(all("/" not in r for r in rows),
      f"row keys must not contain '/' -- MultiQC strips everything before it; got {rows[:2]}")
check(rows[0].startswith("GAPDH |"),
      f"gene symbols must be resolved and the only 2-screen pair first; got {rows[0]!r}")
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
# Candidates renders FIRST in the chimera group (report_section_order), so the
# guide it points at is below it. "... evidence above" sent readers the wrong
# way on a real report.
check("weigh this evidence</strong> below" in desc,
      "Candidates must point DOWN to 'How to weigh this evidence' -- it renders first")
check(all("guide above" not in h.get("description", "")
          for h in d.get("headers", {}).values()),
      "no column tooltip may say 'the guide above' -- the guide is below the table")
for banned in ("ranked by", "confidence tier", "highest confidence"):
    check(banned not in desc.lower(),
          f"section must not present itself as a ranking ({banned!r})")
# nothing here is truncated (5 pairs, nowhere near either threshold), so the
# alphabetical-accident wording must not appear.
check("alphabetically first" not in desc,
      f"no truncation happened here -- 'alphabetically first' must not appear: {desc!r}")
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

# --- Row selection: a real 4-sample run had 66 pairs tied at Screens = 3
# against the OLD fixed cap of 50 -- the report silently dropped 16 of them
# purely alphabetically (genes P-Z). The new rule: show every pair tied at
# the top Screens value; if that's smaller than a useful table, fill down
# whole next-Screens-value tiers; a hard --max-rows cap only ever bites
# when the top group ALONE is bigger than it -- the one case where which
# rows are shown is an alphabetical accident, and the only case where the
# section may say "alphabetically first".
mkdir -p "$T/sel"
cols="gene_id\tte_id\tte_subfamily\tte_family\tte_class\tfound_by\tn_screens\tscreen_evidence\tn_screen_evidence\tcorroboration\tn_corroboration\tcr_events\tcr_reads\tcr_max_samples\tcr_canonical\tcr_chimera_types\ttelocal_active\ttelocal_count\ttelocal_locus\tassembly_transcripts\tassembly_chimera_types\tassembly_strand_match\tassembly_transcript_ids\tsj_events\tsj_reads\tsj_max_samples\tsj_canonical\tsj_chimera_types"

python3 - "$T/sel" <<'PY' || FAIL=1
import gzip, sys
T = sys.argv[1]
cols = ["gene_id","te_id","te_subfamily","te_family","te_class","found_by","n_screens",
        "screen_evidence","n_screen_evidence","corroboration","n_corroboration",
        "cr_events","cr_reads","cr_max_samples","cr_canonical","cr_chimera_types",
        "telocal_active","telocal_count","telocal_locus","assembly_transcripts",
        "assembly_chimera_types","assembly_strand_match","assembly_transcript_ids",
        "sj_events","sj_reads","sj_max_samples","sj_canonical","sj_chimera_types"]

def row(gid, tid, screens):
    d = dict.fromkeys(cols, ".")
    d["gene_id"] = gid; d["te_id"] = tid; d["n_screens"] = str(screens)
    d["found_by"] = "cr"; d["cr_events"] = "1"; d["cr_reads"] = "1"
    d["cr_max_samples"] = "1"; d["cr_canonical"] = "no"
    return "\t".join(d[c] for c in cols)

def write(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write(r + "\n")

# Scenario A: top group (60, screens=3) is bigger than the OLD 50-row cap
# but under --max-rows -- all 60 must be shown, nothing truncated.
write(f"{T}/A.tsv.gz",
      [row(f"G{i:03d}", "T1", 3) for i in range(60)] +
      [row(f"H{i:03d}", "T1", 1) for i in range(5)])

# Scenario B: top group (3, screens=3) is smaller than MIN_ROWS -- must
# fill down with the WHOLE next tier (60 pairs at screens=1).
write(f"{T}/B.tsv.gz",
      [row(f"G{i:03d}", "T1", 3) for i in range(3)] +
      [row(f"H{i:03d}", "T1", 1) for i in range(60)])

# Scenario C: top group (15, screens=3) alone exceeds a small --max-rows
# (10) -- only that case may cut alphabetically into the top group.
write(f"{T}/C.tsv.gz", [row(f"G{i:03d}", "T1", 3) for i in range(15)])
PY

if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/sel/A.tsv.gz" --max-rows 500 \
      --out "$T/sel/A.json" > "$T/sel/A.log" 2>&1; then
  echo "ERROR: scenario A failed"; cat "$T/sel/A.log"; FAIL=1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/sel/B.tsv.gz" --max-rows 500 \
      --out "$T/sel/B.json" > "$T/sel/B.log" 2>&1; then
  echo "ERROR: scenario B failed"; cat "$T/sel/B.log"; FAIL=1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/sel/C.tsv.gz" --max-rows 10 \
      --out "$T/sel/C.json" > "$T/sel/C.log" 2>&1; then
  echo "ERROR: scenario C failed"; cat "$T/sel/C.log"; FAIL=1
fi

python3 - "$T/sel" <<'PY' || FAIL=1
import json, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False

# --- Scenario A: all 60 of the top group shown, nothing truncated -------
a = json.load(open(f"{T}/A.json"))
check(len(a["data"]) == 60, f"scenario A: expected 60 rows shown, got {len(a['data'])}")
desc_a = a.get("description", "")
check("All <strong>60</strong>" in desc_a or "All 60" in desc_a,
      f"scenario A: description must plainly say all 60 pairs are shown: {desc_a!r}")
check("alphabetically first" not in desc_a,
      f"scenario A: nothing truncated, must not say 'alphabetically first': {desc_a!r}")

# --- Scenario B: top group (3) fills down with the whole next tier (60) -
b = json.load(open(f"{T}/B.json"))
check(len(b["data"]) == 63,
      f"scenario B: expected 3 + 60 = 63 rows shown (fill-down), got {len(b['data'])}")
desc_b = b.get("description", "")
check(("All <strong>3</strong>" in desc_b or "All 3" in desc_b)
      and ("60" in desc_b),
      f"scenario B: description must disclose the 3-pair top group and the "
      f"60-pair fill-down with numbers: {desc_b!r}")
check("alphabetically first" not in desc_b,
      f"scenario B: fill-down adds a WHOLE tier, not a partial one -- must "
      f"not say 'alphabetically first': {desc_b!r}")

# --- Scenario C: top group (15) truncated to --max-rows (10) ------------
c = json.load(open(f"{T}/C.json"))
check(len(c["data"]) == 10,
      f"scenario C: expected exactly 10 rows (the --max-rows cap), got {len(c['data'])}")
desc_c = c.get("description", "")
check("15 pairs share the top Screens value (3)" in desc_c,
      f"scenario C: description must name the true top-group size (15) and "
      f"its Screens value (3): {desc_c!r}")
check("alphabetically first 10" in desc_c,
      f"scenario C: description must disclose the cap (10) as the alphabetical "
      f"cutoff: {desc_c!r}")
# the 10 shown must be the alphabetically-first 10 of the 15 (G000..G009)
expected_c = sorted(f"G{i:03d} | T1" for i in range(15))[:10]
check(list(c["data"]) == expected_c,
      f"scenario C: shown rows must be the alphabetically first 10; got "
      f"{list(c['data'])}, want {expected_c}")

sys.exit(0 if ok else 1)
PY

exit $FAIL
