#!/usr/bin/env bash
# Guard 18: empty GTF warning
#
# Run on its own:   .tests/guards/18_empty_gtf_warning.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- an empty-but-existing GTF must produce a parse-time warning
# (not an error) during dry-run.
touch "$T/empty.gtf"
sed -e "s|  gtf: .*|  gtf: $T/empty.gtf|" \
    config/test.yaml > "$T/empty_gtf.yaml"
if ! snakemake --configfile "$T/empty_gtf.yaml" -n --cores 1 >"$T/empty_gtf.log" 2>&1; then
  echo "ERROR: empty GTF dry-run failed"; tail -20 "$T/empty_gtf.log"; FAIL=1
elif ! grep -q "no feature lines" "$T/empty_gtf.log"; then
  echo "ERROR: empty GTF did not produce expected warning"; tail -20 "$T/empty_gtf.log"; FAIL=1
fi

exit $FAIL
