#!/usr/bin/env bash
# Guard 28: classify_chimera_assembly.py classification logic
#
# Run on its own:   .tests/guards/28_classify_chimera_assembly_py_classification_logic.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- synthetic assembly covering all 5 output classes plus a
# negative control (TE on the last exon, no earlier gene match --
# nothing to terminate, must be dropped, not just low-confidence).
printf "chr1\t1000\t1200\tGENE1\t.\t+\n" > "$T/genes.bed"
printf "chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\n" > "$T/exons.bed"
printf "chr1\t500\t700\tTE_A\t.\t+\tL1\tLINE\tL1PA2\nchr1\t2600\t2800\tTE_B\t.\t+\tAluY\tSINE\tAluYa5\nchr1\t9000\t9100\tTE_C\t.\t+\tL2\tLINE\tL2a\nchr1\t1500\t1600\tTE_D\t.\t+\tERV1\tLTR\tMER41\nchr1\t8700\t8800\tTE_E\t.\t+\tL1\tLINE\tL1MA4\n" > "$T/te.bed"
printf 'chr1\tStringTie\ttranscript\t501\t1200\t.\t+\t.\ttranscript_id "T1"; gene_id "MSTRG.1";\n' > "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t501\t700\t.\t+\t.\ttranscript_id "T1"; gene_id "MSTRG.1";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T1"; gene_id "MSTRG.1";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\ttranscript\t1001\t2800\t.\t+\t.\ttranscript_id "T2"; gene_id "MSTRG.2";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T2"; gene_id "MSTRG.2";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t2601\t2800\t.\t+\t.\ttranscript_id "T2"; gene_id "MSTRG.2";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\ttranscript\t1001\t2200\t.\t+\t.\ttranscript_id "T3"; gene_id "MSTRG.3";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T3"; gene_id "MSTRG.3";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t1501\t1600\t.\t+\t.\ttranscript_id "T3"; gene_id "MSTRG.3";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t2001\t2200\t.\t+\t.\ttranscript_id "T3"; gene_id "MSTRG.3";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\ttranscript\t9001\t9100\t.\t+\t.\ttranscript_id "T4"; gene_id "MSTRG.4";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t9001\t9100\t.\t+\t.\ttranscript_id "T4"; gene_id "MSTRG.4";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\ttranscript\t8000\t8800\t.\t+\t.\ttranscript_id "T5"; gene_id "MSTRG.5";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t8000\t8100\t.\t+\t.\ttranscript_id "T5"; gene_id "MSTRG.5";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t8700\t8800\t.\t+\t.\ttranscript_id "T5"; gene_id "MSTRG.5";\n' >> "$T/stringtie.gtf"
if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" --te "$T/te.bed" \
      --breakpoint-tolerance 0 --out "$T/candidates.tsv" > "$T/classify.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/classify.log"; FAIL=1
else
  declare -A want=( [T1]=te_initiated [T2]=te_terminated [T3]=te_exonized [T4]=unspliced_te_only )
  for tid in "${!want[@]}"; do
    got=$(awk -F'\t' -v id="$tid" '$1==id{print $NF}' "$T/candidates.tsv")
    if [ "$got" != "${want[$tid]}" ]; then
      echo "ERROR: $tid classified as '$got', expected '${want[$tid]}'"
      cat "$T/candidates.tsv"; FAIL=1
    fi
  done
  if grep -qP "^T5\t" "$T/candidates.tsv"; then
    echo "ERROR: T5 (TE on last exon, no earlier gene match) should be dropped, not reported"
    cat "$T/candidates.tsv"; FAIL=1
  fi
  te_exon=$(awk -F'\t' '$1=="T2"{print $8"-"$9}' "$T/candidates.tsv")
  if [ "$te_exon" != "2600-2800" ]; then
    echo "ERROR: T2's te_exon_start/end should be 2600-2800 (the TE-overlapping LAST exon, not the first), got $te_exon"
    FAIL=1
  fi
fi

exit $FAIL
