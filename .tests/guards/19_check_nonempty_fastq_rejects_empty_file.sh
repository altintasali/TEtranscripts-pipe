#!/usr/bin/env bash
# Guard 19: check_nonempty_fastq rejects empty file
#
# Run on its own:   .tests/guards/19_check_nonempty_fastq_rejects_empty_file.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- check_nonempty_fastq.py must exit 1 for an empty FASTQ
echo -n "" | gzip > "$T/empty.fq.gz"
if python3 workflow/scripts/check_nonempty_fastq.py "$T/empty.fq.gz" >/dev/null 2>&1; then
  echo "ERROR: check_nonempty_fastq accepted an empty FASTQ"; FAIL=1
fi

exit $FAIL
