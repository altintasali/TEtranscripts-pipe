#!/usr/bin/env bash
# Guard 32: raw FastQC gets a canonical sample name even for oddly-named single-lane fastqs
#
# Run on its own:   .tests/guards/32_raw_fastqc_gets_a_canonical_sample_name_even_for_oddly_nam.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# A single-lane sample's raw fastq keeps its own (possibly
# arbitrary) basename -- merged_fastq_path() skips the merge step
# to avoid an unnecessary copy for the common case. FastQC bakes
# the INPUT's own basename into its report content (the
# "Filename" line in fastqc_data.txt), which is what MultiQC's
# fastqc module uses for the sample name, not the renamed output
# .zip's filename. This checks fastqc_raw's fix (symlink to a
# canonical name before running fastqc) actually lands in the
# report content, not just the DAG -- a dry-run wouldn't catch
# this, since the bug was in report content, not the output path.
mkdir -p "$T/weird_reads"
ln -sf "$(pwd)/.tests/reads/control_rep2_R1.fastq.gz" "$T/weird_reads/SampleA_S1_L001_R1_001.fastq.gz"
{ head -1 .tests/samples.csv; echo "guardsample,$T/weird_reads/SampleA_S1_L001_R1_001.fastq.gz,,auto,control"; } > "$T/weird_samples.csv"
if ! snakemake --configfile config/test.yaml --config samples="$T/weird_samples.csv" --cores 1 \
      results/fastqc/raw/guardsample_R1_fastqc.zip >"$T/weird_fastqc.log" 2>&1; then
  echo "ERROR: fastqc_raw failed for the oddly-named-fastq sample"; tail -40 "$T/weird_fastqc.log"; FAIL=1
else
  FN=$(unzip -p results/fastqc/raw/guardsample_R1_fastqc.zip "*/fastqc_data.txt" | awk -F'\t' '$1=="Filename"{print $2}')
  if [ "$FN" != "guardsample_R1.fastq.gz" ]; then
    echo "ERROR: FastQC's embedded Filename is '$FN', not the canonical guardsample_R1.fastq.gz -- MultiQC General Stats would mismatch"
    FAIL=1
  fi
fi

exit $FAIL
