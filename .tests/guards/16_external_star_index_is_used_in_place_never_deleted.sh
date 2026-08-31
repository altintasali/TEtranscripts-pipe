#!/usr/bin/env bash
# Guard 16: external STAR index is used in place, never deleted
#
# Run on its own:   .tests/guards/16_external_star_index_is_used_in_place_never_deleted.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- The core safety property of build_index: false -- the user's
# external index is a plain INPUT.  No rule may declare it as an
# output (Snakemake wipes a directory() output before the job runs)
# and cleanup_star_index must not touch it even with
# keep_star_index: false.  Fake index: STAR loads SA/SAindex/Genome,
# and those are what common.smk checks for.
fixture_fake_star_index
sed -e 's|^  index: .*|  index: '"$T"'/fake_star_index|' \
    config/test.yaml > "$T/ksif.yaml"
sed -e '/^  index: /a\  build_index: false' "$T/ksif.yaml" > "$T/ksif2.yaml"
mv "$T/ksif2.yaml" "$T/ksif.yaml"
sed -e '/^outputs:/a\'$'\n''  keep_star_index: false' "$T/ksif.yaml" > "$T/ksif2.yaml"
mv "$T/ksif2.yaml" "$T/ksif.yaml"
if ! snakemake --configfile "$T/ksif.yaml" -n --cores 1 -p \
    >"$T/ksif.log" 2>&1; then
  echo "ERROR: keep_star_index=false dry-run failed"; tail -20 "$T/ksif.log"; FAIL=1
elif grep -q "^rule star_index:" "$T/ksif.log"; then
  echo "ERROR: star_index rule planned although build_index=false"; FAIL=1
elif grep -qE "rm -rf '?[^ ']*fake_star_index" "$T/ksif.log"; then
  # Deliberately anchored to the token right after rm -rf: a loose
  # "rm -rf.*fake_star_index" also matches STAR's own
  # "rm -rf <outTmpDir> ... --genomeDir <index>" command line.
  echo "ERROR: external STAR index reached an rm -rf"; FAIL=1
elif ! grep -q "fake_star_index" "$T/ksif.log"; then
  echo "ERROR: external STAR index not used by the alignment jobs"; FAIL=1
fi

# --- ...and with build_index: true the index we built ourselves IS
# cleaned up, so the knob still does what it says.
sed -e '/^outputs:/a\'$'\n''  keep_star_index: false' \
    config/test.yaml > "$T/ksit.yaml"
if ! snakemake --configfile "$T/ksit.yaml" -n --cores 1 -p \
    >"$T/ksit.log" 2>&1; then
  echo "ERROR: build_index=true keep_star_index=false dry-run failed"
  tail -20 "$T/ksit.log"; FAIL=1
elif ! grep -q "cleanup_star_index" "$T/ksit.log"; then
  echo "ERROR: cleanup_star_index rule not planned when we built the index"; FAIL=1
fi

exit $FAIL
