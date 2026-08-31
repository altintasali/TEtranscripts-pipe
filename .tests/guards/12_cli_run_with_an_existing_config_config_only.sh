#!/usr/bin/env bash
# Guard 12: CLI run with an existing config (config-only)
#
# Run on its own:   .tests/guards/12_cli_run_with_an_existing_config_config_only.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- CLI run with an existing --config: no other flags needed -------
# The sheet is taken from the config's 'samples:' key
# (config/test.yaml -> .tests/samples.csv), resolved repo-root
# relative like snakemake itself does.
python workflow/scripts/tetranscripts-pipe run \
  --config config/test.yaml --cores 2 --dry-run > "$T/cli_config_only.log" 2>&1
CLI_RC=$?
if [ $CLI_RC -ne 0 ]; then
  echo "ERROR: CLI run --config-only dry-run failed (exit $CLI_RC)"; tail -40 "$T/cli_config_only.log"; FAIL=1
elif ! grep -q "dry-run OK" "$T/cli_config_only.log"; then
  echo "ERROR: CLI run --config-only did not report success"; tail -40 "$T/cli_config_only.log"; FAIL=1
elif ! grep -q "Job stats:" "$T/cli_config_only.log" || ! grep -q "tecount_summary" "$T/cli_config_only.log"; then
  echo "ERROR: CLI run --config-only planned no jobs (tecount_summary missing)"; tail -40 "$T/cli_config_only.log"; FAIL=1
fi

exit $FAIL
