#!/usr/bin/env bash
# Guard 1: missing reference
#
# Run on its own:   .tests/guards/01_missing_reference.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- a missing/typo'd reference is a hard parse-time error ---------
sed 's|  fasta: .*|  fasta: /nonexistent/ref.fa|' config/test.yaml > "$T/missing_ref.yaml"
if snakemake --configfile "$T/missing_ref.yaml" -n --cores 1 >"$T/missing_ref.log" 2>&1; then
  echo "ERROR: missing reference did not fail the run"; FAIL=1
elif ! grep -q "Missing reference file(s)" "$T/missing_ref.log"; then
  echo "ERROR: wrong failure message for missing reference"; FAIL=1
fi

exit $FAIL
