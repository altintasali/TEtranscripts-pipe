#!/usr/bin/env bash
# Guard 15: build_index=false: index must look like an index
#
# Run on its own:   .tests/guards/15_build_index_false_index_must_look_like_an_index.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- build_index=false pointing at a directory that exists but has
# no STAR index in it is rejected at parse time, not hours later
# when the first alignment job fails inside STAR.
mkdir -p "$T/empty_index"
sed -e 's|^  index: .*|  index: '"$T"'/empty_index\n  build_index: false|' \
    config/test.yaml > "$T/notanindex.yaml"
if snakemake --configfile "$T/notanindex.yaml" -n --cores 1 \
    >"$T/notanindex.log" 2>&1; then
  echo "ERROR: build_index=false with an empty index dir did not fail"; FAIL=1
elif ! grep -q "does not look like a STAR index" "$T/notanindex.log"; then
  echo "ERROR: wrong failure message for non-index directory"; FAIL=1
fi

exit $FAIL
