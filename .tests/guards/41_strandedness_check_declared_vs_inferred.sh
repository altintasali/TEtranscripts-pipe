#!/usr/bin/env bash
# Guard 41: strandedness check: declared vs inferred
#
# Run on its own:   .tests/guards/41_strandedness_check_declared_vs_inferred.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- strandedness.check_provided (default true) runs RSeQC on
# samples that DECLARED a value too, so a sample sheet that
# disagrees with the data gets flagged.  The declared value must
# still win for quantification -- reporting a mismatch is the
# point; silently overriding the user's sheet is not.
mkdir -p "$T/strand/qc"
awk -F, 'BEGIN{OFS=","} NR==1{print;next}
  {if($1=="control_rep2") $4="auto"; else $4="forward"; print}' \
  .tests/samples.csv > "$T/strand/mixed.csv"
sed "s|^samples: .*|samples: $T/strand/mixed.csv|" config/test.yaml \
  > "$T/strand/on.yaml"
if ! snakemake --configfile "$T/strand/on.yaml" -n --cores 1 -p \
    >"$T/strand/on.log" 2>&1; then
  echo "ERROR: mixed declared/auto dry-run failed"
  tail -20 "$T/strand/on.log"; FAIL=1
else
  N=$(sed -n '/^Job stats:/,/^$/p' "$T/strand/on.log" \
      | awk '$1=="rseqc_infer_experiment"{print $2; exit}')
  if [ "$N" != "6" ]; then
    echo "ERROR: check_provided should run RSeQC on all 6 samples, got $N"; FAIL=1
  fi
  if ! grep -q "strandedness_check" "$T/strand/on.log"; then
    echo "ERROR: strandedness_check rule not planned"; FAIL=1
  fi
  if ! grep -q -- "--stranded forward" "$T/strand/on.log"; then
    echo "ERROR: declared strandedness did not reach TEcount"; FAIL=1
  fi
fi
# ...and check_provided: false runs RSeQC only where it is needed.
printf '\nstrandedness:\n  min_fraction: 0.6\n  check_provided: false\n' \
  >> "$T/strand/on.yaml"
if ! snakemake --configfile "$T/strand/on.yaml" -n --cores 1 \
    >"$T/strand/off.log" 2>&1; then
  echo "ERROR: check_provided=false dry-run failed"
  tail -20 "$T/strand/off.log"; FAIL=1
else
  N=$(sed -n '/^Job stats:/,/^$/p' "$T/strand/off.log" \
      | awk '$1=="rseqc_infer_experiment"{print $2; exit}')
  if [ "$N" != "1" ]; then
    echo "ERROR: check_provided=false should run RSeQC once (the auto sample), got $N"; FAIL=1
  fi
fi

# --- the table itself: OK / MISMATCH / auto-detected, and it must
# render in a real MultiQC run next to the built-in RSeQC plot.
mkstrand() {
  printf '\n\nThis is PairEnd Data\nFraction of reads failed to determine: %s\nFraction of reads explained by "1++,1--,2+-,2-+": %s\nFraction of reads explained by "1+-,1-+,2++,2--": %s\n' \
    "$2" "$3" "$4" > "$T/strand/qc/$1_infer_experiment.txt"
  printf '%s\n' "$5" > "$T/strand/$1_strandedness.txt"
}
mkstrand okS   0.05 0.03 0.92 reverse
mkstrand badS  0.05 0.03 0.92 reverse
mkstrand autoS 0.04 0.91 0.05 forward
python3 - "$T/strand" <<'PY' || FAIL=1
# run_name must be __main__: these scripts guard their main() call so
# the unit tests can import them, and snakemake's script: directive
# executes them as __main__ too.
import builtins, json, runpy, sys, types
class NS(dict):
    def __getattr__(self, k): return self[k]
T = sys.argv[1]
smp = ["okS", "badS", "autoS"]
out = f"{T}/qc/strandedness_check_mqc.json"
builtins.snakemake = types.SimpleNamespace(
    params=NS(samples=smp,
              declared={"okS": "reverse", "badS": "forward", "autoS": "auto"}),
    input=NS(reports=[f"{T}/qc/{s}_infer_experiment.txt" for s in smp],
             calls=[f"{T}/{s}_strandedness.txt" for s in smp]),
    output=[out])
runpy.run_path("workflow/scripts/strandedness_check_mqc.py", run_name="__main__")
d = json.load(open(out))
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
check("data" in d, "custom_content needs a top-level data key")
rows = d["data"]
check(rows["okS"]["status"] == "OK", "agreeing sample must be OK")
check(rows["badS"]["status"] == "MISMATCH", "disagreeing sample must be MISMATCH")
check(rows["badS"]["used"] == "forward", "MISMATCH must still USE the declared value")
check(rows["autoS"]["status"] == "auto-detected", "auto sample mislabelled")
check(rows["autoS"]["used"] == "forward", "auto sample must use the inferred value")
sys.exit(0 if ok else 1)
PY
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/strand/out" -n r "$T/strand/qc" > "$T/strand/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the strandedness check"
  tail -30 "$T/strand/render.log"; FAIL=1
elif ! grep -q "Strandedness check" "$T/strand/out/r.html"; then
  echo "ERROR: strandedness check table not rendered"; FAIL=1
elif ! grep -q "MISMATCH" "$T/strand/out/r.html"; then
  echo "ERROR: mismatch not surfaced in the report"; FAIL=1
elif ! grep -q "Infer experiment" "$T/strand/out/r.html"; then
  echo "ERROR: built-in RSeQC infer_experiment plot missing"; FAIL=1
fi

exit $FAIL
