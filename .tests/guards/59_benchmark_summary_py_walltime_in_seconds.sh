#!/usr/bin/env bash
# Guard 59: benchmark_summary.py reports wall time as plain seconds
# (mean/max, two numeric columns), and the "Resource Usage" section
# always renders (no config gate -- reverted after 2026 feedback that it
# should just be there, not opt-in).
#
# BUG FIXED 2026, TWICE: this column was first six fixed h/min/s mean/max
# columns, replaced with one adaptive string column ("45.2s (max
# 3.20min)", unit picked per value) because a fixed hour column rounded a
# fast rule to "0.000". That fix introduced a second bug: a STRING column
# sorts lexicographically, not by duration, so clicking the header to find
# the slowest rule silently did the wrong thing, and comparing two rules
# meant reading past mismatched units by eye -- exactly what a user
# reported. Plain seconds (not the old six h/min/s columns) fixes both:
# one column pair, always comparable, never rounds a fast rule away.
#
# Run on its own:   .tests/guards/59_benchmark_summary_py_walltime_in_seconds.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- part 1: script-level output shape --------------------------------
mkdir -p "$T/bm/star_index/star_index"
# Fast rule (sub-second) and a slow one (multi-minute) -- both must render
# in the SAME unit (seconds), not rounded away or unit-switched.
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
    # the rule names its benchmark list (input.benchmarks) since it also
    # declares its own script as an input -- see guard 65
    input=types.SimpleNamespace(benchmarks=[sys.argv[1], sys.argv[2]],
                                script="workflow/scripts/benchmark_summary.py"),
    params=NS(allocated={"star_index": {"threads": 4, "mem_mb": 8000},
                          "star_align": {"threads": 2, "mem_mb": 16000}}),
    output=[out],
)
runpy.run_path("workflow/scripts/benchmark_summary.py", run_name="__main__")
doc = json.load(open(out))
ok = True
headers = doc["headers"]
# ANTI-REGRESSION: neither the old six-fixed-column shape nor the
# intermediate one-adaptive-string-column shape may come back.
old_six = ["walltime_mean_h", "walltime_max_h", "walltime_mean_min",
           "walltime_max_min"]
leaked = [c for c in old_six if c in headers]
if leaked:
    print(f"ERROR: old fixed-unit (min/h) walltime columns present: {leaked}"); ok = False
if "walltime" in headers:
    print("ERROR: old single adaptive-string 'walltime' column still present"); ok = False
for want in ("walltime_mean_s", "walltime_max_s"):
    if want not in headers:
        print(f"ERROR: {want!r} column missing from headers"); ok = False

fast = doc["data"].get("star_index", {})
slow = doc["data"].get("star_align", {})
# both plain numbers, in SECONDS, not strings and not unit-switched
if fast.get("walltime_mean_s") != 0.8 or fast.get("walltime_max_s") != 0.8:
    print(f"ERROR: fast rule's walltime wrong: {fast}"); ok = False
if slow.get("walltime_mean_s") != 245.0 or slow.get("walltime_max_s") != 245.0:
    print(f"ERROR: slow rule's walltime wrong (must stay in seconds, not minutes): {slow}"); ok = False
for row in (fast, slow):
    for key in ("walltime_mean_s", "walltime_max_s"):
        if isinstance(row.get(key), str):
            print(f"ERROR: {key} is a string ({row[key]!r}), not numeric -- "
                  "a string column sorts lexicographically, not by duration")
            ok = False
sys.exit(0 if ok else 1)
PY
if [ "${FAIL:-0}" != "0" ]; then
  echo "ERROR: benchmark_summary.py walltime-in-seconds check failed"
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
