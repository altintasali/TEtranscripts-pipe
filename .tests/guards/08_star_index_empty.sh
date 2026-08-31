#!/usr/bin/env bash
# Guard 8: star.index empty
#
# Run on its own:   .tests/guards/08_star_index_empty.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- star.index empty string: same default applies --------------------
sed 's|^  index: .*|  index: ""|' config/test.yaml > "$T/empty_index.yaml"
if ! snakemake --configfile "$T/empty_index.yaml" -n --cores 1 >"$T/empty_index.log" 2>&1; then
  echo "ERROR: empty star.index must dry-run (default results/star_index)"; FAIL=1
fi

exit $FAIL
