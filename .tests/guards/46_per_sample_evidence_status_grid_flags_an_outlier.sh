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
# An HTML table, not a heatmap. The heatmap plotted the status ENCODING
# (1.0/0.5/0.0), so a healthy cohort rendered as a wall of "1" -- a number
# meaning nothing to the reader, with the counts behind the verdict hidden.
check(d["plot_type"] == "html", "status grid must be an html table")
body = d["data"]
import re
text = re.sub(r"<[^>]+>", " ", body)
text = re.sub(r"\s+", " ", text)

for s in ["GV_KO_01", "GV_KO_02", "GV_WT_01", "GV_WT_02"]:
    check(s in text, f"row for {s} missing")
check("Strandedness" in text, "strandedness column missing")

# the counts themselves must be on the page -- that is the whole point
check("118,000" in text, "cohort-scale count not shown (GV_KO_01 events_total)")
check("11,800" in text, "the 10x-low sample's own count not shown")
check("cohort median" in text, "median reference row missing")

# the outlier must be flagged red, its peers green
def cells(sample):
    m = re.search(r"<tr><td><strong>" + sample + r"</strong>(.*?)</tr>", body, re.S)
    check(m is not None, f"could not locate the {sample} row")
    return re.findall(r"background:(#[0-9a-f]{6})", m.group(1)) if m else []

low = cells("GV_WT_01")
check(low.count("#fbe6e6") >= 8, f"the 10x-low sample must be red across yield layers, got {low}")
for good in ("GV_KO_01", "GV_KO_02", "GV_WT_02"):
    g = cells(good)
    check(all(c == "#e4f3e4" for c in g), f"{good} is within 2x of the median and must be green, got {g}")

check("cohort median" in d["description"], "description must say the scale is relative")
check("fall outside 2x" in text, "the verdict line must state how many cells are flagged")
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
