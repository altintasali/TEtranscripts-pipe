#!/usr/bin/env bash
# Guard 26: chimera.assembly.enabled: false plans nothing
#
# Run on its own:   .tests/guards/26_chimera_assembly_enabled_false_plans_nothing.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- the opt-out direction: turning the screen off must remove
# every one of its rules from the DAG, not merely skip them.
sed '/^  assembly:$/{n;n;s/enabled: true/enabled: false/}' \
    config/test.yaml > "$T/assembly_off.yaml"
if ! grep -qE "^    enabled: false" "$T/assembly_off.yaml"; then
  echo "ERROR: guard 26 sed did not flip assembly.enabled to false"
  sed -n '/^  assembly:/,+3p' "$T/assembly_off.yaml"; FAIL=1
elif ! snakemake --configfile "$T/assembly_off.yaml" -n --cores 1 >"$T/assembly_off.log" 2>&1; then
  echo "ERROR: chimera.assembly.enabled: false must dry-run"; tail -40 "$T/assembly_off.log"; FAIL=1
elif grep -qE "stringtie_assemble|chimera_assembly_classify|star_align_for_assembly" "$T/assembly_off.log"; then
  echo "ERROR: chimera_assembly rules still planned with enabled: false"
  tail -40 "$T/assembly_off.log"; FAIL=1
fi

exit $FAIL
