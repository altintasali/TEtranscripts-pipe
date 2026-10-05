#!/usr/bin/env bash
# Guard 57: the renamed chimera.sj_junctions config key is rejected with a fix
#
# Run on its own:   .tests/guards/57_chimera_sj_junctions_config_key_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- chimera.sj_junctions was renamed to chimera.splice_junctions:
# "sj_junctions" repeated itself ("SJ" already means splice junction), and
# the config-key abbreviation was hard to read even though the pipeline
# already used the short "sj"/"SJ" form pervasively elsewhere. Same
# reasoning/pattern as guards 52/56, and it sits right next to them in
# common/config.smk.
sed 's/^  splice_junctions:$/  sj_junctions:/' config/test.yaml > "$T/old_key.yaml"
if ! grep -q "^  sj_junctions:" "$T/old_key.yaml"; then
  echo "ERROR: fixture did not produce an old-style chimera.sj_junctions key"; FAIL=1
elif snakemake --configfile "$T/old_key.yaml" -n --cores 1 > "$T/old.log" 2>&1; then
  echo "ERROR: a config using chimera.sj_junctions must NOT parse"; FAIL=1
else
  for probe in "chimera.sj_junctions" "chimera.splice_junctions"; do
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
  echo "ERROR: chimera.splice_junctions (the current key) must parse"; tail -20 "$T/new.log"; FAIL=1
fi

exit $FAIL
