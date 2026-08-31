#!/usr/bin/env bash
# Guard 13: qc-off keeps the summary barplots
#
# Run on its own:   .tests/guards/13_qc_off_keeps_the_summary_barplots.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- tetranscripts.qc.enabled: false still plans the summary barplots
# The assignment + TE-class barplots use the raw cntTables, so they
# are decoupled from the R-based sample-QC view (v0.6.5). The 4-space
# indent pins this to tetranscripts.qc.enabled (not trimming/chimera).
sed 's/^    enabled: true/    enabled: false/' config/test.yaml > "$T/qc_off.yaml"
if ! snakemake --configfile "$T/qc_off.yaml" -n --cores 1 >"$T/qc_off.log" 2>&1; then
  echo "ERROR: qc-off config must dry-run"; tail -40 "$T/qc_off.log"; FAIL=1
elif ! grep -q "tecount_summary" "$T/qc_off.log"; then
  echo "ERROR: tecount_summary not planned with tetranscripts.qc.enabled false"; tail -40 "$T/qc_off.log"; FAIL=1
elif grep -qE "tecount_counts|tecount_qc_transform" "$T/qc_off.log"; then
  echo "ERROR: R-based tecount QC rules still planned with tetranscripts.qc.enabled false"; tail -40 "$T/qc_off.log"; FAIL=1
fi

exit $FAIL
