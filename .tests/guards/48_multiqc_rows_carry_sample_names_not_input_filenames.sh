#!/usr/bin/env bash
# Guard 48: MultiQC rows carry sample names, not input filenames
#
# Run on its own:   .tests/guards/48_multiqc_rows_carry_sample_names_not_input_filenames.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- On a real run General Statistics grew four extra rows named
# after the sequencing facility's FASTQ files, plus a separate
# {sample}_flagstat row per sample. Two different causes:
#
#   flagstat  -- filename suffix, fixed by extra_fn_clean_exts.
#   trimming  -- TrimGalore! writes "Input filename: ..." INSIDE the
#                report and MultiQC's cutadapt module reads the
#                sample name from that line, so renaming the report
#                does nothing. The input is symlinked to a canonical
#                {sample}_R{n} name and TrimGalore! runs on that.
#
# Single-lane samples are handed their RAW path, so only they were
# affected -- the bundled test data hid it. This fixture is
# deliberately facility-named and single-lane.
mkdir -p "$T/nm/raw" "$T/nm/trimming" "$T/nm/samtools"
cp .tests/reads/treatment_rep3_R1.fastq.gz "$T/nm/raw/MUX11129_01_WT_GV1.fastq.gz"
( set -eu; cd "$T/nm"
  stage=trimming/.stage_GV_WT_01; rm -rf "$stage"; mkdir -p "$stage"
  set --; i=1
  # shellcheck disable=SC2043  # one read on purpose: this mirrors the
  # pipeline's per-read loop for a single-end sample
  for f in raw/MUX11129_01_WT_GV1.fastq.gz; do
    case "$f" in *.gz) ext=.fastq.gz;; *) ext=.fastq;; esac
    ln -s "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" \
          "$stage/GV_WT_01_R${i}$ext"
    set -- "$@" "$stage/GV_WT_01_R${i}$ext"; i=$((i+1))
  done
  trim_galore --gzip --cores 1 --fastqc_args "-t 1" --basename GV_WT_01 \
    --output_dir trimming "$@" > tg.log 2>&1
  rm -rf "$stage" ) || { echo "ERROR: staged trim_galore run failed"; FAIL=1; }
printf '39900000 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 secondary\n0 + 0 supplementary\n0 + 0 duplicates\n39900000 + 0 mapped (100.00%% : N/A)\n0 + 0 paired in sequencing\n' \
  > "$T/nm/samtools/GV_WT_01_flagstat.txt"
mkdir -p "$T/nm/mq"
cp "$T/nm/trimming"/*_trimming_report.txt "$T/nm/trimming"/*_fastqc.zip \
   "$T/nm/samtools"/*_flagstat.txt "$T/nm/mq/" 2>/dev/null || true
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/nm/out" -n r "$T/nm/mq" > "$T/nm/mq.log" 2>&1; then
  echo "ERROR: multiqc failed on the naming fixture"; tail -20 "$T/nm/mq.log"; FAIL=1
else
  python3 - "$T/nm/out/r_data/multiqc_data.json" <<'PY' || FAIL=1
import json, sys
d = json.load(open(sys.argv[1]))
names = {n for v in d.get("report_data_sources", {}).values()
           for s in v.values() for n in s}
ok = True
if any("MUX11129" in n for n in names):
    print(f"ERROR: an input filename leaked into MultiQC sample names: {sorted(names)}")
    ok = False
if any(n.endswith("_flagstat") for n in names):
    print(f"ERROR: _flagstat suffix not cleaned: {sorted(names)}")
    ok = False
if "GV_WT_01" not in names:
    print(f"ERROR: expected a GV_WT_01 row, got {sorted(names)}")
    ok = False
sys.exit(0 if ok else 1)
PY
fi

exit $FAIL
