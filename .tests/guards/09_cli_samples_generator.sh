#!/usr/bin/env bash
# Guard 9: CLI samples generator
#
# Run on its own:   .tests/guards/09_cli_samples_generator.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- the samples CLI must reproduce the committed fixture -----------
# (same rows, order-insensitive; generated paths are absolute while
# the fixture's are relative, so compare basenames)
python workflow/scripts/tetranscripts-pipe samples --reads .tests/reads --dry-run > "$T/cli_sheet.csv" 2> "$T/cli_samples.err"
if [ $? -ne 0 ]; then
  echo "ERROR: samples generator failed"; cat "$T/cli_samples.err"; FAIL=1
else
  norm() { grep -v '^#' "$1" | grep -v '^sample,' | cut -d, -f1-3 | sed -E 's|,[^,]*/|,|g' | sort; }
  norm "$T/cli_sheet.csv" > "$T/cli_norm.txt"
  norm .tests/samples.csv > "$T/fixture_norm.txt"
  if ! diff "$T/cli_norm.txt" "$T/fixture_norm.txt" > /dev/null; then
    echo "ERROR: samples generator output differs from the committed fixture"
    diff "$T/cli_norm.txt" "$T/fixture_norm.txt"
    FAIL=1
  fi
fi

exit $FAIL
