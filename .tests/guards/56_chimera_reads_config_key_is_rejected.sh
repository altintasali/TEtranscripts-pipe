#!/usr/bin/env bash
# Guard 56: the renamed chimera.reads config key is rejected with a fix
#
# Run on its own:   .tests/guards/56_chimera_reads_config_key_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- chimera.reads was renamed to chimera.chimeric_reads: "reads" was
# ambiguous (the SJ.out.tab screen is also read-derived), and this screen's
# report label, config key and data columns used to be three different
# words for the same thing. Same reasoning/pattern as guard 52's
# chimera.junction check, and it sits right next to it in common/config.smk.
sed 's/^  chimeric_reads:$/  reads:/' config/test.yaml > "$T/old_key.yaml"
if ! grep -q "^  reads:" "$T/old_key.yaml"; then
  echo "ERROR: fixture did not produce an old-style chimera.reads key"; FAIL=1
elif snakemake --configfile "$T/old_key.yaml" -n --cores 1 > "$T/old.log" 2>&1; then
  echo "ERROR: a config using chimera.reads must NOT parse"; FAIL=1
else
  for probe in "chimera.reads" "chimera.chimeric_reads"; do
    if ! grep -q -- "$probe" "$T/old.log"; then
      echo "ERROR: the error message never mentions '$probe'"; FAIL=1
    fi
  done
  if grep -q "ValidationError" "$T/old.log" && ! grep -q "was renamed to" "$T/old.log"; then
    echo "ERROR: fell through to the bare jsonschema error instead of the migration message"
    FAIL=1
  fi
fi

# ...and the current key still works, so the check cannot fire spuriously.
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/new.log" 2>&1; then
  echo "ERROR: chimera.chimeric_reads (the current key) must parse"; tail -20 "$T/new.log"; FAIL=1
fi

exit $FAIL
