#!/usr/bin/env bash
# Guard 7: star.index omitted
#
# Run on its own:   .tests/guards/07_star_index_omitted.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- star.index omitted: defaults to generated results/star_index -----
sed '/^  index: .*star_index$/d' config/test.yaml > "$T/no_index.yaml"
if ! snakemake --configfile "$T/no_index.yaml" -n --cores 1 >"$T/no_index.log" 2>&1; then
  echo "ERROR: omitted star.index must dry-run (default results/star_index)"; FAIL=1
elif ! grep -q "results/star_index" "$T/no_index.log"; then
  echo "ERROR: default star index path results/star_index not used"; FAIL=1
fi

exit $FAIL
