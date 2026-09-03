#!/usr/bin/env bash
# Guard 52: the renamed chimera.junction config key is rejected with a fix
#
# Run on its own:   .tests/guards/52_chimera_junction_config_key_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- chimera.junction was renamed to chimera.reads in 0.12.0. The schema is
# additionalProperties: false, so the old key fails on its own -- but with a
# bare jsonschema error naming "junction" and nothing else. Every user with an
# existing input/config.yaml hits that on upgrade, so common/config.smk checks
# for it FIRST and says what to rename. This pins the message, not just the
# failure: a schema error alone would satisfy "it fails" while leaving the
# user stuck.
sed 's/^  reads:$/  junction:/' config/test.yaml > "$T/old_key.yaml"
if ! grep -q "^  junction:" "$T/old_key.yaml"; then
  echo "ERROR: fixture did not produce an old-style chimera.junction key"; FAIL=1
elif snakemake --configfile "$T/old_key.yaml" -n --cores 1 > "$T/old.log" 2>&1; then
  echo "ERROR: a config using chimera.junction must NOT parse"; FAIL=1
else
  for probe in "chimera.junction" "chimera.reads" "0.12.0"; do
    if ! grep -q -- "$probe" "$T/old.log"; then
      echo "ERROR: the error message never mentions '$probe'"; FAIL=1
    fi
  done
  # it must be OUR actionable error, not the raw schema rejection
  if grep -q "ValidationError" "$T/old.log" && ! grep -q "was renamed to" "$T/old.log"; then
    echo "ERROR: fell through to the bare jsonschema error instead of the migration message"
    FAIL=1
  fi
fi

# ...and the current key still works, so the check cannot fire spuriously.
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/new.log" 2>&1; then
  echo "ERROR: chimera.reads (the current key) must parse"; tail -20 "$T/new.log"; FAIL=1
fi

exit $FAIL
