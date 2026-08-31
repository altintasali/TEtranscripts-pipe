#!/usr/bin/env bash
# Guard 46: per-sample evidence status grid flags an outlier
#
# Run on its own:   .tests/guards/46_per_sample_evidence_status_grid_flags_an_outlier.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The cohort-level heatmaps cannot answer "is ONE sample
# different"; a replicate at a tenth of its peers' yield is absorbed
# by any summary. Fixture makes GV_WT_01 exactly that, and the grid
# must flag it while leaving the other three green. Thresholds are
# relative to the cohort median, so this also pins that absolute
# cut-offs are not being used.
mkdir -p "$T/st/qc"
python3 - "$T/st" <<'PY'
import gzip, json, sys
T = sys.argv[1]
def qc(sample, scale):
    m = {"events_total": 118000*scale, "events_gene_te": 5300*scale,
         "canonical_gene_to_te": 530*scale, "canonical_te_to_gene": 220*scale,
         "ambiguous_direction": 400*scale,
         "chimera_type_te_initiated": 900*scale,
         "chimera_type_te_exonized": 3500*scale,
         "chimera_type_te_terminated": 700*scale,
         "strand_match_yes": 2600*scale}
    with gzip.open(f"{T}/{sample}_junction_qc.tsv.gz", "wt") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"sample\t{sample}\n")
        for k, v in m.items():
            fh.write(f"{k}\t{int(v)}\n")
for s, sc in [("GV_KO_01",1.0),("GV_KO_02",1.1),("GV_WT_01",0.1),("GV_WT_02",0.95)]:
    qc(s, sc)
json.dump({"data": {s: {"status": "OK"} for s in
                    ["GV_KO_01","GV_KO_02","GV_WT_01","GV_WT_02"]}},
          open(f"{T}/strand.json", "w"))
PY
if ! python3 workflow/scripts/sample_evidence_status_mqc.py \
    --qc-tables "$T/st"/GV_KO_01_junction_qc.tsv.gz \
                "$T/st"/GV_KO_02_junction_qc.tsv.gz \
                "$T/st"/GV_WT_01_junction_qc.tsv.gz \
                "$T/st"/GV_WT_02_junction_qc.tsv.gz \
    --samples GV_KO_01 GV_KO_02 GV_WT_01 GV_WT_02 \
    --strandedness "$T/st/strand.json" \
    --out "$T/st/qc/sample_evidence_status_mqc.json" \
    > "$T/st/log" 2>&1; then
  echo "ERROR: sample_evidence_status_mqc.py failed"; cat "$T/st/log"; FAIL=1
else
  python3 - "$T/st" <<'PY' || FAIL=1
import json, sys
d = json.load(open(f"{sys.argv[1]}/qc/sample_evidence_status_mqc.json"))
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
check(d["plot_type"] == "heatmap", "status grid must be a heatmap")
check(d["ycats"] == ["GV_KO_01","GV_KO_02","GV_WT_01","GV_WT_02"], "rows must be samples")
check("Strandedness" in d["xcats"], "strandedness column missing")
check(len(d["data"]) == 4 and all(len(r) == len(d["xcats"]) for r in d["data"]),
      "grid shape does not match its labels")
check(all(v in (0.0, 0.5, 1.0) for r in d["data"] for v in r), "cells must be pass/warn/fail")
i = d["ycats"].index("GV_WT_01")
yield_cols = [c for c, n in enumerate(d["xcats"]) if n != "Strandedness"]
check(all(d["data"][i][c] == 0.0 for c in yield_cols),
      "the 10x-low sample must fail every yield layer")
for good in ("GV_KO_01", "GV_KO_02", "GV_WT_02"):
    g = d["ycats"].index(good)
    check(all(d["data"][g][c] == 1.0 for c in yield_cols),
          f"{good} is within 2x of the median and must pass")
check("cohort median" in d["description"], "description must say the scale is relative")
sys.exit(0 if ok else 1)
PY
fi
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/st/out" -n r "$T/st/qc" > "$T/st/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the status grid"
  tail -30 "$T/st/render.log"; FAIL=1
elif ! grep -q "Per-sample status" "$T/st/out/r.html"; then
  echo "ERROR: per-sample status grid not rendered"; FAIL=1
fi

exit $FAIL
