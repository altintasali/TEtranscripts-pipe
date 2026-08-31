#!/usr/bin/env bash
# Guard 20: telocal enabled with empty locind auto-builds
#
# Run on its own:   .tests/guards/20_telocal_enabled_with_empty_locind_auto_builds.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- telocal.enabled: true with empty locind must now succeed
# (auto-builds from TE GTF via telocal_locind rule)
printf '\ntelocal:\n  enabled: true\n  locind: ""\n' >> "$T/telocal_empty_locind.yaml"
cat config/test.yaml "$T/telocal_empty_locind.yaml" > "$T/telocal_auto.yaml"
if ! snakemake --configfile "$T/telocal_auto.yaml" -n --cores 1 >"$T/telocal_auto.log" 2>&1; then
  echo "ERROR: telocal enabled with empty locind failed (auto-build should work)"
  tail -20 "$T/telocal_auto.log"
  FAIL=1
elif ! grep -q "telocal_locind" "$T/telocal_auto.log"; then
  echo "ERROR: telocal_locind rule not triggered in auto-build mode"
  FAIL=1
fi

exit $FAIL
