#!/usr/bin/env bash
# Guard 47: MultiQC ignores unparseable bulk but keeps its content
#
# Run on its own:   .tests/guards/47_multiqc_ignores_unparseable_bulk_but_keeps_its_content.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- MultiQC's search directories are results/* subtrees that also
# hold the pipeline's data: trimmed FASTQs, BAMs, gzipped counts
# matrices, every rule log. On a real run the search phase stalled
# long enough to look like a hang. fn_ignore_* now skips them --
# this pins that the skips do not also drop the custom content
# MultiQC is there to render.
mkdir -p "$T/ign/fastqc/raw" "$T/ign/telocal/qc" \
         "$T/ign/pipeline_info/logs" "$T/ign/pipeline_info/benchmarks/star_align" \
         "$T/ign/trimming"
gzip -dc .tests/reads/treatment_rep2_R1.fastq.gz > "$T/ign/a.fastq"
fastqc -o "$T/ign/fastqc/raw" -q "$T/ign/a.fastq" > /dev/null 2>&1
printf '{"id":"keeper","section_name":"Keeper","plot_type":"html","data":"<p>kept</p>"}\n' \
  > "$T/ign/pipeline_info/keeper_mqc.json"
printf '{"id":"tlkeeper","section_name":"TElocal Keeper","plot_type":"html","data":"<p>kept</p>"}\n' \
  > "$T/ign/telocal/qc/telocal_keeper_mqc.json"
# decoy: valid custom content INSIDE an ignored path -- if
# fn_ignore_paths works it must never be rendered.
printf '{"id":"decoy","section_name":"DecoySection","plot_type":"html","data":"<p>no</p>"}\n' \
  > "$T/ign/pipeline_info/logs/decoy_mqc.json"
head -c 3000000 /dev/urandom > "$T/ign/trimming/s_trimmed.fq.gz"
awk 'BEGIN{printf "feature\ta\tb\n"; for(i=0;i<200000;i++) printf "L1_dup%d\t1\t2\n", i}' \
  | gzip -c > "$T/ign/telocal/qc/log2_counts.tsv.gz"
printf 'log line\n' > "$T/ign/pipeline_info/logs/rule.log"
printf 's\th:m:s\nx\t0:00:01\n' > "$T/ign/pipeline_info/benchmarks/star_align/s.txt"
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
    -o "$T/ign/out" -n r "$T/ign" > "$T/ign/log" 2>&1; then
  echo "ERROR: multiqc failed with the ignore rules"; tail -30 "$T/ign/log"; FAIL=1
elif ! grep -q "Keeper" "$T/ign/out/r.html"; then
  echo "ERROR: custom content in pipeline_info/ was dropped"; FAIL=1
elif ! grep -q "TElocal Keeper" "$T/ign/out/r.html"; then
  echo "ERROR: custom content beside an ignored matrix was dropped"; FAIL=1
elif grep -q "DecoySection" "$T/ign/out/r.html"; then
  echo "ERROR: fn_ignore_paths did not exclude pipeline_info/logs"; FAIL=1
elif ! grep -q "FastQC" "$T/ign/out/r.html"; then
  echo "ERROR: fastqc reports were dropped"; FAIL=1
fi

exit $FAIL
