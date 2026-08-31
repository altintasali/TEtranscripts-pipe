#!/usr/bin/env bash
# Guard 24: star.two_pass default (cohort) pools junctions across samples
#
# Run on its own:   .tests/guards/24_star_two_pass_default_cohort_pools_junctions_across_sample.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- default config/test.yaml (no override -- star.two_pass
# defaults to "cohort") must plan star_align_pass1 for every
# sample, star_merge_junctions once, and hand the merged file to
# star_align via --sjdbFileChrStartEnd.
if ! snakemake --configfile config/test.yaml -n -p --cores 1 --until star_align >"$T/two_pass_cohort.log" 2>&1; then
  echo "ERROR: default (cohort) star.two_pass must dry-run"; tail -40 "$T/two_pass_cohort.log"; FAIL=1
elif ! grep -q "star_align_pass1" "$T/two_pass_cohort.log"; then
  echo "ERROR: star_align_pass1 not planned by default (star.two_pass should default to cohort)"
  tail -40 "$T/two_pass_cohort.log"; FAIL=1
elif ! grep -q "star_merge_junctions" "$T/two_pass_cohort.log"; then
  echo "ERROR: star_merge_junctions not planned by default"
  tail -40 "$T/two_pass_cohort.log"; FAIL=1
elif ! grep -q -- "--sjdbFileChrStartEnd results/star_pass1/merged_SJ.out.tab" "$T/two_pass_cohort.log"; then
  echo "ERROR: star_align missing --sjdbFileChrStartEnd with default (cohort) star.two_pass"
  tail -40 "$T/two_pass_cohort.log"; FAIL=1
fi

# --- merge_splice_junctions.py: pools reads across samples so a
# junction too weak in any single sample survives if the combined
# support clears the threshold, drops annotated junctions and ones
# that never clear it even pooled.
printf "chr1\t1000\t2000\t1\t2\t0\t1\t0\t20\nchr1\t5000\t6000\t1\t2\t1\t10\t0\t20\n" > "$T/mgA_SJ.out.tab"
printf "chr1\t1000\t2000\t1\t2\t0\t1\t0\t20\n" > "$T/mgB_SJ.out.tab"
printf "chr1\t9000\t9500\t2\t1\t0\t1\t0\t15\n" > "$T/mgC_SJ.out.tab"
if ! python3 workflow/scripts/merge_splice_junctions.py \
      --sj-tables "$T/mgA_SJ.out.tab" "$T/mgB_SJ.out.tab" "$T/mgC_SJ.out.tab" \
      --min-unique-reads 2 --out "$T/merged_SJ.out.tab" >"$T/merge_sj.log" 2>&1; then
  echo "ERROR: merge_splice_junctions.py failed"; cat "$T/merge_sj.log"; FAIL=1
elif [ "$(cat "$T/merged_SJ.out.tab")" != "$(printf 'chr1\t1000\t2000\t1\t2\t0\t2\t0\t20')" ]; then
  echo "ERROR: merge_splice_junctions.py did not pool/filter junctions as expected"
  echo "-- got --"; cat "$T/merged_SJ.out.tab"; FAIL=1
fi

exit $FAIL
