#!/usr/bin/env bash
# Guard 29: chimera_assembly_summary_mqc.py renders in MultiQC
#
# Run on its own:   .tests/guards/29_chimera_assembly_summary_mqc_py_renders_in_multiqc.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

printf "transcript_id\tgtf_gene_id\tchrom\tstrand\tn_exons\ttranscript_start\ttranscript_end\tte_exon_start\tte_exon_end\tte_id\tte_family\tte_class\tte_overlap_exon_rank\tte_hits_all\tmatched_gene_id\tmatched_gene_strand\tgene_hits_all\tstrand_match\tchimera_type\tconfirmed_by_junction_screen\tjunction_supporting_reads\n" > "$T/ca.tsv"
printf "T1\tMSTRG.1\tchr1\t+\t2\t500\t1200\t500\t700\tTE_A\tL1\tLINE\t1\tTE_A\tGENE1\t+\tGENE1\tyes\tte_initiated\tyes\t7\n" >> "$T/ca.tsv"
gzip -c "$T/ca.tsv" > "$T/ca.tsv.gz"
mkdir -p "$T/custom/chimera/qc"
if ! python workflow/scripts/chimera_assembly_summary_mqc.py \
      --candidates "$T/ca.tsv.gz" \
      --out-classes "$T/custom/chimera/qc/chimera_assembly_classes_mqc.json" \
      --out-highlights "$T/custom/chimera/qc/chimera_assembly_highlights_mqc.json" \
      > "$T/casum.log" 2>&1; then
  echo "ERROR: chimera_assembly_summary_mqc.py failed"; cat "$T/casum.log"; FAIL=1
fi
if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/mqc_ca" -n report "$T/custom" > "$T/mqc_ca.log" 2>&1; then
  echo "ERROR: multiqc run with chimera_assembly failed"; tail -40 "$T/mqc_ca.log"; FAIL=1
else
  for cid in chimera_assembly_classes chimera_assembly_highlights; do
    if ! grep -q "custom_content | $cid: Found" "$T/mqc_ca.log"; then
      echo "ERROR: section '$cid' not rendered by MultiQC"
      grep "custom_content" "$T/mqc_ca.log" | head -20
      FAIL=1
    fi
  done
  if ! grep -q "What it cannot see" "$T/mqc_ca/report.html"; then
    echo "ERROR: assembly blind-spot note missing from rendered report"
    FAIL=1
  fi
fi

exit $FAIL
