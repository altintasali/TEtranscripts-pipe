#!/usr/bin/env bash
# Guard 39: 'TE analysis' evidence overview reflects what ran
#
# Run on its own:   .tests/guards/39_start_here_evidence_overview_reflects_what_ran.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The report's reading guide is generated from the resolved
# switches, so its central claim -- that the two chimera screens
# are the only independent pair -- must track them.
mkdir -p "$T/eo/qc"
python3 - "$T/eo/qc/evidence_overview_mqc.json" <<'PY' || FAIL=1
# run_name must be __main__: these scripts guard their main() call so
# the unit tests can import them, and snakemake's script: directive
# executes them as __main__ too.
import builtins, json, runpy, sys, types
class NS(dict):
    def __getattr__(self, k): return self[k]
out = sys.argv[1]
ok = True
cases = [
    (dict(_telocal_enabled=True, _chimera_junction_enabled=True,
          _chimera_assembly_enabled=True, _two_pass="cohort",
          _has_condition=True), "Two independent"),
    (dict(_telocal_enabled=False, _chimera_junction_enabled=True,
          _chimera_assembly_enabled=False, _two_pass="none",
          _has_condition=False), "One chimera screen"),
    (dict(_telocal_enabled=False, _chimera_junction_enabled=False,
          _chimera_assembly_enabled=False, _two_pass="per_sample",
          _has_condition=False), "quantification only"),
]
for params, expected in cases:
    builtins.snakemake = types.SimpleNamespace(
        params=NS(_sample_count=6, **params), output=[out])
    runpy.run_path("workflow/scripts/evidence_overview_mqc.py", run_name="__main__")
    doc = json.load(open(out))
    if "data" not in doc:
        print("ERROR: custom_content needs a top-level data key"); ok = False
    if expected not in doc["data"]:
        print(f"ERROR: expected {expected!r} for {params}"); ok = False
sys.exit(0 if ok else 1)
PY
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/eo/out" -n r "$T/eo/qc" > "$T/eo/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the evidence overview section"
  tail -30 "$T/eo/render.log"; FAIL=1
elif ! grep -q "TE analysis" "$T/eo/out/r.html"; then
  echo "ERROR: 'TE analysis' section not rendered"; FAIL=1
fi

exit $FAIL
