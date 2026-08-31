#!/usr/bin/env bash
# Guard 6: gz refs + empty decompressed_dir
#
# Run on its own:   .tests/guards/06_gz_refs_empty_decompressed_dir.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_gz_refs

# --- gz refs + empty decompressed_dir: hard parse error ---------------
sed 's|  sjdb_overhang: 49|  sjdb_overhang: 49\n  decompressed_dir: ""|' "$T/gz.yaml" > "$T/gz_empty.yaml"
if snakemake --configfile "$T/gz_empty.yaml" -n --cores 1 >"$T/gz_empty.log" 2>&1; then
  echo "ERROR: gz refs + empty decompressed_dir did not fail the run"; FAIL=1
elif ! grep -q "(empty)" "$T/gz_empty.log"; then
  echo "ERROR: wrong failure message for gz refs + empty decompressed_dir"; FAIL=1
fi

exit $FAIL
