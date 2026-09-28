#!/usr/bin/env bash
# Guard 59: benchmark_summary.py reports ONE adaptive wall-time column
# (not six fixed h/min/s mean/max columns), and the "Resource Usage"
# section always renders (no config gate -- reverted after 2026 feedback
# that it should just be there, not opt-in).
#
# Run on its own:   .tests/guards/59_benchmark_summary_py_adaptive_walltime_column.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- part 1: script-level output shape --------------------------------
mkdir -p "$T/bm/star_index/star_index"
# Fast rule (sub-second) and a slow one (multi-minute), so the adaptive
# formatter must pick a DIFFERENT unit for each, in the same table.
printf "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n" > "$T/bm/star_index/star_index.txt"
printf "0.8\t0:00:00\t120\t150\t100\t110\t0\t0\t95.0\t0.76\n" >> "$T/bm/star_index/star_index.txt"
mkdir -p "$T/bm/star_align"
printf "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n" > "$T/bm/star_align/S1.txt"
printf "245.0\t0:04:05\t8192\t9000\t8000\t8100\t0\t0\t380.0\t930.0\n" >> "$T/bm/star_align/S1.txt"

python3 - "$T/bm/star_index/star_index.txt" "$T/bm/star_align/S1.txt" "$T/out.json" <<'PY' || FAIL=1
import builtins, json, runpy, sys, types
class NS(dict):
    def __getattr__(self, k): return self[k]
out = sys.argv[3]
builtins.snakemake = types.SimpleNamespace(
    input=[sys.argv[1], sys.argv[2]],
    params=NS(allocated={"star_index": {"threads": 4, "mem_mb": 8000},
                          "star_align": {"threads": 2, "mem_mb": 16000}}),
    output=[out],
)
runpy.run_path("workflow/scripts/benchmark_summary.py", run_name="__main__")
doc = json.load(open(out))
ok = True
headers = doc["headers"]
old_cols = ["walltime_mean_h", "walltime_max_h", "walltime_mean_min",
            "walltime_max_min", "walltime_mean_s", "walltime_max_s"]
leaked = [c for c in old_cols if c in headers]
if leaked:
    print(f"ERROR: old fixed-unit walltime columns still present: {leaked}"); ok = False
if "walltime" not in headers:
    print("ERROR: single adaptive 'walltime' column missing from headers"); ok = False
fast = doc["data"].get("star_index", {}).get("walltime", "")
slow = doc["data"].get("star_align", {}).get("walltime", "")
if "0.8s" not in fast:
    print(f"ERROR: sub-second rule should render in seconds, got {fast!r}"); ok = False
if "min" not in slow:
    print(f"ERROR: multi-minute rule should render in minutes, got {slow!r}"); ok = False
sys.exit(0 if ok else 1)
PY
if [ "${FAIL:-0}" != "0" ]; then
  echo "ERROR: benchmark_summary.py adaptive-walltime check failed"
fi

# --- part 2: Resource Usage always renders, no config key needed ------
# ANTI-REGRESSION: this section was briefly gated behind
# outputs.report_resource_usage (default off) -- reverted; benchmark_summary
# must be unconditionally in the planned DAG, and the config key must be
# gone from the schema (a stale key nobody reads is worse than no key).
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/plan.log" 2>&1; then
  echo "ERROR: dry run failed"; tail -30 "$T/plan.log"; FAIL=1
elif ! grep -qE "^benchmark_summary\b" "$T/plan.log"; then
  echo "ERROR: benchmark_summary rule not scheduled -- Resource Usage must always render"
  FAIL=1
fi

python3 - <<'PY' || FAIL=1
import yaml
schema = yaml.safe_load(open("workflow/schemas/config.schema.yaml"))
prop = schema["properties"]["outputs"]["properties"].get("report_resource_usage")
if prop is not None:
    print("ERROR: outputs.report_resource_usage should be gone from config.schema.yaml, found it still declared")
    raise SystemExit(1)
PY

exit $FAIL
