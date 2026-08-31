#!/usr/bin/env bash
# Guard 14: CLI with config in a different checkout (deployed split)
#
# Run on its own:   .tests/guards/14_cli_with_config_in_a_different_checkout_deployed_split.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- CLI installed elsewhere, config in an analysis checkout -------
# The deployment runs this script from a shared apps/ copy while
# pointing --config at an analysis working copy. Relative paths in
# the config (samples: input/samples.csv) must then resolve against
# the analysis checkout, and snakemake must run from that root.
WD2="$T/deployed_wd"
mkdir -p "$WD2"
cp -r workflow config .tests "$WD2/"
mkdir -p "$WD2/input"
sed 's|samples: .*|samples: input/samples.csv|' config/test.yaml > "$WD2/input/config.yaml"
cp .tests/samples.csv "$WD2/input/samples.csv"
python workflow/scripts/tetranscripts-pipe run \
  --config "$WD2/input/config.yaml" --cores 2 --dry-run > "$T/deployed.log" 2>&1
DEPLOY_RC=$?
if [ $DEPLOY_RC -ne 0 ]; then
  echo "ERROR: CLI run with a config in another checkout failed (exit $DEPLOY_RC)"; tail -40 "$T/deployed.log"; FAIL=1
elif ! grep -q "dry-run OK" "$T/deployed.log"; then
  echo "ERROR: deployed-split run did not report success"; tail -40 "$T/deployed.log"; FAIL=1
elif ! grep -q "$WD2" "$T/deployed.log"; then
  echo "ERROR: run did not anchor to the config's checkout ($WD2)"; tail -40 "$T/deployed.log"; FAIL=1
fi

exit $FAIL
