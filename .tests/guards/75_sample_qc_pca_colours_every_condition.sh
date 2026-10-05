#!/usr/bin/env bash
# Guard 75: the sample-QC PCA gives every condition its own colour
#
# sample_qc.R coloured PCA points from a fixed 7-colour palette indexed by
# group number, so groups 8+ got the colour "NA" and looked identical. A
# real 84-sample run with 10 conditions (genotype_stage) showed 8 distinct
# "colours" for 10 groups. It now keeps the Okabe-Ito colours up to 7
# groups and switches to an evenly spaced HCL palette beyond that.
#
# Run on its own:   .tests/guards/75_sample_qc_pca_colours_every_condition.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[guard 75] SKIP: Rscript not on PATH"; exit 0
fi

# 12 samples, 12 conditions (more than the 7 fixed colours), 30 features
python3 - "$T" <<'PY'
import random, sys
T = sys.argv[1]
random.seed(1)
samples = [f"S{i:02d}" for i in range(12)]
with open(f"{T}/mat.tsv", "w") as fh:
    fh.write("feature\t" + "\t".join(samples) + "\n")
    for f in range(30):
        fh.write(f"f{f}\t" + "\t".join(f"{random.gauss(5, 1):.3f}" for _ in samples) + "\n")
with open(f"{T}/samples.csv", "w") as fh:
    fh.write("sample,condition\n")
    for i, s in enumerate(samples):
        fh.write(f"{s},cond{i:02d}\n")
PY

if ! Rscript workflow/scripts/sample_qc.R --plots tecount "$T/mat.tsv" "$T/samples.csv" \
      1 vst "$T/pca.json" "$T/heat.json" > "$T/r.log" 2>&1; then
  echo "ERROR: sample_qc.R --plots failed"; cat "$T/r.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import json, re, sys
T = sys.argv[1]
d = json.load(open(f"{T}/pca.json"))
colours = [v["color"] for v in d["data"].values()]
ok = True
if any(not re.fullmatch(r"#[0-9A-Fa-f]{6}", c) for c in colours):
    print(f"ERROR: every PCA point needs a hex colour (no 'NA'); got {sorted(set(colours))}")
    ok = False
if len(set(colours)) != 12:
    print(f"ERROR: 12 conditions need 12 distinct colours; got {len(set(colours))}")
    ok = False
sys.exit(0 if ok else 1)
PY

exit $FAIL
