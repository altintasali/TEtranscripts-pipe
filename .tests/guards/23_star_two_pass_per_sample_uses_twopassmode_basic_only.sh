#!/usr/bin/env bash
# Guard 23: star.two_pass: per_sample uses --twopassMode Basic only
#
# Run on its own:   .tests/guards/23_star_two_pass_per_sample_uses_twopassmode_basic_only.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- per_sample must add --twopassMode Basic to star_align and add
# no other rules (no pass1/merge stage -- that is cohort-only).
sed '/^  index: .tests\/resources\/star_index$/a\  two_pass: per_sample' \
    config/test.yaml > "$T/two_pass_per_sample.yaml"
if ! snakemake --configfile "$T/two_pass_per_sample.yaml" -n -p --cores 1 --until star_align >"$T/two_pass_per_sample.log" 2>&1; then
  echo "ERROR: star.two_pass: per_sample must dry-run"; tail -40 "$T/two_pass_per_sample.log"; FAIL=1
elif ! grep -q -- "--twopassMode Basic" "$T/two_pass_per_sample.log"; then
  echo "ERROR: --twopassMode Basic missing from star_align with star.two_pass: per_sample"
  tail -40 "$T/two_pass_per_sample.log"; FAIL=1
elif grep -qE "star_align_pass1|star_merge_junctions" "$T/two_pass_per_sample.log"; then
  echo "ERROR: star_align_pass1/star_merge_junctions planned with star.two_pass: per_sample"
  tail -40 "$T/two_pass_per_sample.log"; FAIL=1
fi

exit $FAIL
