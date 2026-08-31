#!/usr/bin/env bash
# Guard 30: telocal.indexer defaults to fast
#
# Run on its own:   .tests/guards/30_telocal_indexer_defaults_to_fast.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- default config/test.yaml (no override) must plan
# build_telocal_index.py with --indexer fast.
if ! snakemake --configfile config/test.yaml -n -p --cores 1 --until telocal_locind >"$T/indexer_default.log" 2>&1; then
  echo "ERROR: default config must dry-run"; tail -40 "$T/indexer_default.log"; FAIL=1
elif ! grep -q -- "--indexer fast" "$T/indexer_default.log"; then
  echo "ERROR: telocal_locind should default to --indexer fast"
  tail -40 "$T/indexer_default.log"; FAIL=1
fi

exit $FAIL
