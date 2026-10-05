#!/usr/bin/env bash
# Guard 66: the report records which commit of the pipeline built it
#
# VERSION ("0.14.2") does not change between commits, so a report built from
# an older checkout -- or one whose report steps never re-ran after a fix --
# looked exactly like a current one. PIPELINE_GIT_STATE (rules/common/
# envs.smk) records the git commit plus whether workflow/ has uncommitted
# changes, and rule config_used passes it as a PARAM so a new commit re-runs
# config_used and MultiQC.
#
# This runs the real _pipeline_git_state() (extracted from envs.smk, not a
# copy) in three throwaway checkouts: clean, dirty, and not-a-git-repo.
#
# Run on its own:   .tests/guards/66_report_records_the_pipeline_commit.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- wiring: config_used takes it as a param, the script shows it ---------
if ! sed -n '/^rule config_used:/,/^rule /p' workflow/rules/qc.smk \
     | grep -q '_pipeline_commit=PIPELINE_GIT_STATE'; then
  echo "ERROR: rule config_used must pass PIPELINE_GIT_STATE as a param (the rerun trigger)"
  FAIL=1
fi
if ! grep -q '"pipeline_commit"' workflow/scripts/config_used_mqc.py; then
  echo "ERROR: config_used_mqc.py no longer emits a pipeline_commit row"
  FAIL=1
fi

# --- behaviour: clean / dirty / not a checkout -----------------------------
if ! command -v git >/dev/null 2>&1; then
  echo "[guard 66] SKIP behaviour checks: git not on PATH"
  exit $FAIL
fi

mkdir -p "$T/clean/workflow" "$T/nogit/workflow"
printf 'rule a:\n    output: "x"\n' > "$T/clean/workflow/Snakefile"
cp "$T/clean/workflow/Snakefile" "$T/nogit/workflow/Snakefile"
(
  cd "$T/clean" && git init -q && git add workflow \
    && git -c user.email=g@x -c user.name=guard commit -q -m init
) || { echo "ERROR: could not set up the throwaway git checkout"; exit 1; }
EXPECT=$(git -C "$T/clean" rev-parse --short=12 HEAD)

python3 - "$T" "$EXPECT" <<'PY' || FAIL=1
import os, re, sys
tmp, expect = sys.argv[1], sys.argv[2]
src = open("workflow/rules/common/envs.smk").read()
m = re.search(r"^def _pipeline_git_state\(\):\n(?:(?:    .*|)\n)+", src, re.M)
if not m:
    print("ERROR: _pipeline_git_state() not found in envs.smk"); sys.exit(1)
ns = {"os": os}
exec(m.group(0), ns)
fn = ns["_pipeline_git_state"]

def state_in(d):
    cwd = os.getcwd()
    os.chdir(d)  # the pipeline runs from the repo root, as here
    try:
        return fn()
    finally:
        os.chdir(cwd)

ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False

clean = state_in(os.path.join(tmp, "clean"))
check(clean == expect, f"clean checkout: expected {expect!r}, got {clean!r}")

with open(os.path.join(tmp, "clean", "workflow", "Snakefile"), "a") as fh:
    fh.write("# local edit\n")
dirty = state_in(os.path.join(tmp, "clean"))
check(dirty == f"{expect} + uncommitted changes in workflow/",
      f"edited workflow/: expected the commit plus a dirty marker, got {dirty!r}")

nogit = state_in(os.path.join(tmp, "nogit"))
check(nogit.startswith("unknown"),
      f"outside a git checkout: expected 'unknown (...)', got {nogit!r}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
