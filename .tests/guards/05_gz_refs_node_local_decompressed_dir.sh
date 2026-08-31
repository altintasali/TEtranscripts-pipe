#!/usr/bin/env bash
# Guard 5: gz refs + node-local decompressed_dir
#
# Run on its own:   .tests/guards/05_gz_refs_node_local_decompressed_dir.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_gz_refs

# --- gz refs + node-local /tmp decompressed_dir: hard parse error ----
sed 's|  sjdb_overhang: 49|  sjdb_overhang: 49\n  decompressed_dir: /tmp/rnaseq-decompressed|' "$T/gz.yaml" > "$T/gz_tmpdir.yaml"
if snakemake --configfile "$T/gz_tmpdir.yaml" -n --cores 1 >"$T/gz_tmpdir.log" 2>&1; then
  echo "ERROR: gz refs + node-local decompressed_dir did not fail the run"; FAIL=1
elif ! grep -q "node-local temp directory" "$T/gz_tmpdir.log"; then
  echo "ERROR: wrong failure message for gz refs + node-local decompressed_dir"; FAIL=1
fi

exit $FAIL
