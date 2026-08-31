#!/usr/bin/env bash
# Guard 40: wildcards are pinned to real sample names
#
# Run on its own:   .tests/guards/40_wildcards_are_pinned_to_real_sample_names.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- Without a {sample} constraint, star_align_pass1's
# results/star_pass1/{sample}_SJ.out.tab also matches
# star_merge_junctions' results/star_pass1/merged_SJ.out.tab.
# Snakemake then tries the per-sample rule with sample="merged",
# the input function raises KeyError, and it silently recovers --
# but only in non-strict DAG mode.  Both properties are checked:
# no recovered stack trace, and a clean run under strict mode.
if ! snakemake --configfile config/test.yaml -n --cores 1 \
    >"$T/wc.log" 2>&1; then
  echo "ERROR: default dry-run failed"; tail -20 "$T/wc.log"; FAIL=1
elif grep -q "KeyError: 'merged'" "$T/wc.log"; then
  echo "ERROR: merged_SJ.out.tab still matches the per-sample rule"; FAIL=1
fi
if ! snakemake --configfile config/test.yaml -n --cores 1 \
    --strict-dag-evaluation cyclic-graph functions periodic-wildcards \
    >"$T/wc_strict.log" 2>&1; then
  echo "ERROR: dry-run fails under --strict-dag-evaluation"
  tail -20 "$T/wc_strict.log"; FAIL=1
fi

exit $FAIL
