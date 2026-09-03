#!/usr/bin/env bash
# Guard 51: count views are counts, not percentages
#
# Run on its own:   .tests/guards/51_count_views_are_counts_not_percentages.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- Two plots pair a COUNTS dataset with a RATE dataset behind one
# switcher. MultiQC derives the axis/tooltip suffix from ylab when ysuffix
# is unset ("%" in ylab -> suffix "%", multiqc/plots/plot.py), and a
# PLOT-level ylab is inherited by every dataset -- so a "% ..." ylab
# silently stamped "%" on the counts view too, rendering raw counts as
# "241%" and "1,898%".
#
# The 6000% axis was that same suffix bug -- stacked COUNTS wearing a "%"
# sign -- not the stacking, so the canonical plot stays stacked on purpose
# (its counts view stacks to the sample's total canonical junctions).
# chimera_assembly_strand_rate_plot is the one that must NOT stack: it draws
# one bar per class, each a rate over its own denominator.
#
# This pins the suffix on every dataset of both plots, plus each plot's
# deliberate stacking choice, so neither is silently flipped again.
mkdir -p "$T/p"

# Junction metric tables -- totals large, canonical a modest fraction, so a
# counts view mislabelled "%" is unmistakable (values far above 100).
python3 - "$T" <<'PY'
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
from junction_qc_mqc import DIRECTIONS
T = sys.argv[1]
for i, s in enumerate(["S1", "S2"]):
    with open(f"{T}/p/{s}.tsv", "w") as fh:
        # junction_qc.py writes a "metric\tvalue" header and load_metrics
        # skips it; without one here the FIRST metric is silently eaten.
        fh.write("metric\tvalue\n")
        for j, d in enumerate(DIRECTIONS):
            tot = 500 + 100 * j + 50 * i
            fh.write(f"direction_{d}\t{tot}\n")
            fh.write(f"canonical_{d}\t{int(tot * (0.03 + 0.01 * j))}\n")
PY

python3 - "$T" <<'PY'
import sys, gzip, csv
T = sys.argv[1]
rows = []
for cls, (n, m) in {"te_exonized": (40, 30), "te_initiated": (25, 5),
                    "te_terminated": (18, 12)}.items():
    rows += [{"chimera_type": cls, "gene_id": f"G{i}", "te_id": f"T{i}",
              "strand_match": "yes" if i < m else "no"} for i in range(n)]
with gzip.open(f"{T}/p/candidates.tsv.gz", "wt") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
PY

python3 workflow/scripts/junction_qc_mqc.py \
  --tables "$T/p/S1.tsv" "$T/p/S2.tsv" --samples S1 S2 \
  --out "$T/p/junction_mqc.json" \
  --out-canonical "$T/p/canonical_mqc.json" > "$T/p/emit.log" 2>&1
python3 workflow/scripts/chimera_assembly_summary_mqc.py \
  --candidates "$T/p/candidates.tsv.gz" \
  --out-classes "$T/p/asm_classes_mqc.json" \
  --out-highlights "$T/p/asm_high_mqc.json" \
  --out-strand-rate "$T/p/asm_strand_mqc.json" >> "$T/p/emit.log" 2>&1

if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/p/out" -n r "$T/p" > "$T/p/mqc.log" 2>&1; then
  echo "ERROR: multiqc failed"; tail -30 "$T/p/mqc.log"; FAIL=1
else
  python3 - "$T/p/out/r_data/multiqc_data.json" <<'PY'
import json, sys
plots = json.load(open(sys.argv[1]))["report_plot_data"]
bad = []
for pid in ("chimera_canonical_rate_plot", "chimera_assembly_strand_rate_plot"):
    if pid not in plots:
        bad.append(f"{pid}: plot missing from the report")
        continue
    o = plots[pid]
    barmode = (o.get("layout") or {}).get("barmode")
    if pid == "chimera_assembly_strand_rate_plot" and barmode != "group":
        bad.append(f"{pid}: barmode is {barmode!r}, not 'group' -- one bar per "
                   "class, each a rate over its own denominator, must not stack")
    if pid == "chimera_canonical_rate_plot" and barmode == "group":
        bad.append(f"{pid}: barmode is 'group'; this plot is stacked on "
                   "purpose (the 6000% axis was the suffix bug, not stacking)")
    for ds in o.get("datasets", []):
        label = ds.get("label", "")
        lay = ds.get("layout") or {}
        sufs = {ax: (lay.get(ax) or {}).get("ticksuffix", "")
                for ax in ("xaxis", "yaxis")}
        carried = {a: s for a, s in sufs.items() if s}
        if label.strip().startswith("%"):
            if "%" not in sufs.values():
                bad.append(f"{pid}/{label!r}: rate view lost its '%' suffix")
        elif "%" in sufs.values():
            bad.append(f"{pid}/{label!r}: COUNTS view carries a '%' suffix "
                       f"({carried}) -- raw counts render as '241%'")
if bad:
    print("ERROR: " + "\nERROR: ".join(bad))
    sys.exit(1)
print("count/rate views carry the right suffix; canonical stacks, "
      "assembly strand-rate does not")
PY
  [ $? -eq 0 ] || FAIL=1
fi

exit $FAIL
