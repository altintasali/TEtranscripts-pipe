#!/usr/bin/env bash
# Guard 17: build_index=false without ref.fasta
#
# Run on its own:   .tests/guards/17_build_index_false_without_ref_fasta.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_fake_star_index

# --- build_index=false: ref.fasta is optional.  Remove fasta from
# the config and confirm the dry-run still passes.
sed -e 's|^  index: .*|  index: '"$T"'/fake_star_index|' \
    -e '/^  fasta: /d' \
    config/test.yaml > "$T/no_fasta.yaml"
sed -e '/^  index: /a\  build_index: false' "$T/no_fasta.yaml" > "$T/no_fasta2.yaml"
mv "$T/no_fasta2.yaml" "$T/no_fasta.yaml"
if ! snakemake --configfile "$T/no_fasta.yaml" -n --cores 1 >"$T/no_fasta.log" 2>&1; then
  echo "ERROR: build_index=false without ref.fasta failed"; tail -20 "$T/no_fasta.log"; FAIL=1
elif grep -q "Missing reference" "$T/no_fasta.log"; then
  echo "ERROR: build_index=false without ref.fasta raised missing reference error"; FAIL=1
fi

exit $FAIL
