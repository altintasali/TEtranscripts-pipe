#!/usr/bin/env bash
# Guard 4: gzipped refs, default decompressed_dir
#
# Run on its own:   .tests/guards/04_gzipped_refs_default_decompressed_dir.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- gzipped references: .gz copies of the test refs (lib.sh) --------
fixture_gz_refs

# --- gz refs + default decompressed_dir: dry-run must succeed --------
if ! snakemake --configfile "$T/gz.yaml" -n --cores 1 >"$T/gz.log" 2>&1; then
  echo "ERROR: gzipped refs with default decompressed_dir must dry-run"; FAIL=1
fi

exit $FAIL
