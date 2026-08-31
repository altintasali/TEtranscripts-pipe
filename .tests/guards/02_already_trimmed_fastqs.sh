#!/usr/bin/env bash
# Guard 2: already-trimmed fastqs
#
# Run on its own:   .tests/guards/02_already_trimmed_fastqs.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- already-trimmed-looking fastqs are rejected while trimming on --
printf 'sample,fastq_1,fastq_2,strandedness,condition\nzyg_wt_01,zyg_wt_01_trimmed.fq.gz,,auto,wt\n' > "$T/trimmed.csv"
sed "s|samples: .*|samples: $T/trimmed.csv|" config/test.yaml > "$T/trimmed.yaml"
if snakemake --configfile "$T/trimmed.yaml" -n --cores 1 >"$T/trimmed.log" 2>&1; then
  echo "ERROR: already-trimmed fastqs did not fail the run"; FAIL=1
elif ! grep -q "already look trimmed" "$T/trimmed.log"; then
  echo "ERROR: wrong failure message for already-trimmed fastqs"; FAIL=1
fi

exit $FAIL
