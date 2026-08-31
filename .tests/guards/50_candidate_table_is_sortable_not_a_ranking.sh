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
# and its order is a COUNT the reader can re-sort rather than a score the
# pipeline asserts.
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

# every signal is its own column, so any of them can be sorted on -- Gene and
# TE insertion included, so a reader can find a specific gene
headers = d.get("headers", {})
for col in ("Gene", "TE insertion", "Evidence types", "Splice motif",
            "Samples", "Found by", "Strand match", "Reads"):
    check(col in headers, f"column {col!r} missing -- the reader cannot sort on it")
# the table opens on the evidence count; sort_rows alone does not stick
check(d["pconfig"].get("defaultsort"),
      "defaultsort must set the opening order -- sort_rows does not survive")

# rows keep candidates.tsv.gz's own order (n_evidence desc, then alphabetical)
rows = list(d["data"])
# " | ", never " / ": MultiQC cleans table row names like filenames and splits
# on "/", so a "GENE / te_id" key rendered as just the TE id and the gene
# vanished from the report entirely.
check(all("/" not in r for r in rows),
      f"row keys must not contain '/' -- MultiQC strips everything before it; got {rows[:2]}")
check(rows[0].startswith("GAPDH |"),
      f"gene symbols must be resolved and the densest-evidence pair first; got {rows[0]!r}")
n_ev = [d["data"][r]["Evidence types"] for r in rows]
check(n_ev == sorted(n_ev, reverse=True),
      f"default order must follow n_evidence descending; got {n_ev}")

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

# ...and it must render, with the table toolbox present.
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
sys.exit(0 if ok else 1)
PY2
then
  echo "ERROR: gene names did not survive MultiQC rendering"; FAIL=1
fi

exit $FAIL
