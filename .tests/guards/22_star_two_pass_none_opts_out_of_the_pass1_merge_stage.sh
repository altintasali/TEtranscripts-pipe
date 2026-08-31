#!/usr/bin/env bash
# Guard 22: star.two_pass: none opts out of the pass1/merge stage
#
# Run on its own:   .tests/guards/22_star_two_pass_none_opts_out_of_the_pass1_merge_stage.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- star.two_pass: none must plan the original, single-pass-only
# DAG (star.two_pass defaults to "cohort" -- see common.smk -- so
# this checks that explicitly opting out still works).
sed '/^  index: .tests\/resources\/star_index$/a\  two_pass: none' \
    config/test.yaml > "$T/two_pass_none.yaml"
if ! snakemake --configfile "$T/two_pass_none.yaml" -n --cores 1 >"$T/two_pass_none.log" 2>&1; then
  echo "ERROR: star.two_pass: none must dry-run"; tail -40 "$T/two_pass_none.log"; FAIL=1
elif grep -qE "star_align_pass1|star_merge_junctions" "$T/two_pass_none.log"; then
  echo "ERROR: star_align_pass1/star_merge_junctions planned with star.two_pass: none"
  tail -40 "$T/two_pass_none.log"; FAIL=1
fi

exit $FAIL
