#!/usr/bin/env bash
# Guard 80: rule inputs and outputs do not embed the run directory's path
#
# Snakemake records every job's input paths for its provenance check. The
# workflow's own scripts, the MultiQC config and an in-run STAR index used to
# be declared by ABSOLUTE path, which embeds the run directory's name: on a
# real 84-sample run renamed after it finished, snakemake wanted to rerun
# 1,811 of 2,064 jobs ("set of input files has changed"). They are now
# relative to the run directory, which is where snakemake always runs.
#
# Checked on a dry-run of config/test.yaml from the repo root (the run
# directory here): no input: or output: line may contain the repo root's
# absolute path. Paths a reader sees in the report stay absolute -- the
# config table resolves them (see .tests/unit/test_config_used_mqc.py).
#
# Run on its own:   .tests/guards/80_rule_paths_do_not_embed_the_run_directory.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/dry.log" 2>&1; then
  echo "ERROR: dry-run of config/test.yaml failed"; tail -30 "$T/dry.log"; exit 1
fi
root="$(pwd -P)"
# the dry-run prints each job's input:/output: lines
grep -E '^\s+(input|output): ' "$T/dry.log" > "$T/io.txt"
if [ ! -s "$T/io.txt" ]; then
  echo "ERROR: no input:/output: lines in the dry-run output -- cannot check"; FAIL=1
fi
if grep -qF "$root/" "$T/io.txt" || grep -qF "$GUARD_REPO_ROOT/" "$T/io.txt"; then
  echo "ERROR: rule inputs/outputs embed the run directory's absolute path, so"
  echo "renaming or moving a finished run would make snakemake rerun it:"
  grep -oE "($root|$GUARD_REPO_ROOT)/[^, ]+" "$T/io.txt" | sort -u | head -10
  FAIL=1
fi
# the script declarations are actually there (the check above is not vacuous)
grep -qE 'workflow/scripts/[A-Za-z0-9_]+\.(py|R)' "$T/io.txt" \
  || { echo "ERROR: no workflow/scripts/ inputs found in the dry-run"; FAIL=1; }

exit $FAIL
