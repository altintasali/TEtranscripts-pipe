#!/usr/bin/env bash
# Guard 10: CLI run --dry-run
#
# Run on its own:   .tests/guards/10_cli_run_dry_run.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- CLI run must scaffold config + sheet and pass a real dry-run ---
python workflow/scripts/tetranscripts-pipe run \
  --config "$T/cli/config.yaml" --samplesheet "$T/cli/samples.csv" \
  --reads .tests/reads \
  --fasta .tests/resources/genome.fa \
  --gtf .tests/resources/genome.gtf \
  --te-gtf .tests/resources/te_annotation.gtf \
  --cores 2 --dry-run > "$T/cli_run.log" 2>&1
CLI_RC=$?
if [ $CLI_RC -ne 0 ]; then
  echo "ERROR: CLI run --dry-run failed (exit $CLI_RC)"; tail -40 "$T/cli_run.log"; FAIL=1
elif ! grep -q "dry-run OK" "$T/cli_run.log"; then
  echo "ERROR: CLI run --dry-run did not report success"; tail -40 "$T/cli_run.log"; FAIL=1
elif ! grep -q "Job stats:" "$T/cli_run.log" || ! grep -q "tecount" "$T/cli_run.log"; then
  echo "ERROR: CLI run --dry-run planned no jobs (tecount missing from the DAG)"; tail -40 "$T/cli_run.log"; FAIL=1
fi

exit $FAIL
