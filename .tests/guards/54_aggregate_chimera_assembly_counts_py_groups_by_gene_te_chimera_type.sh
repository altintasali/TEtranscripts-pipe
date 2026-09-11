#!/usr/bin/env bash
# Guard 54: aggregate_chimera_assembly_counts.py groups by
# (matched_gene_id, te_id, chimera_type), not by transcript_id or by
# (gene, te) alone, and drops non-candidate rows.
#
# Run on its own:   .tests/guards/54_aggregate_chimera_assembly_counts_py_groups_by_gene_te_chimera_type.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- transcripts.tsv.gz-shaped fixture (classify_chimera_assembly.py's own
# header). Only transcript_id(1), te_id(10), matched_gene_id(16), and
# chimera_type(20) matter here; the rest are cheap placeholders.
printf "transcript_id\tgtf_gene_id\tchrom\tstrand\tn_exons\ttranscript_start\ttranscript_end\tte_exon_start\tte_exon_end\tte_id\tte_subfamily\tte_family\tte_class\tte_overlap_exon_rank\tte_hits_all\tmatched_gene_id\tmatched_gene_strand\tgene_hits_all\tstrand_match\tchimera_type\n" > "$T/transcripts.tsv"
# T1, T2: same (GENE1, TE_A, te_initiated) -- must be summed together.
printf "T1\tMSTRG.1\tchr1\t+\t2\t1\t100\t1\t50\tTE_A\tL1PA2\tL1\tLINE\t1\tTE_A\tGENE1\t+\tGENE1\tyes\tte_initiated\n" >> "$T/transcripts.tsv"
printf "T2\tMSTRG.2\tchr1\t+\t2\t1\t100\t1\t50\tTE_A\tL1PA2\tL1\tLINE\t1\tTE_A\tGENE1\t+\tGENE1\tyes\tte_initiated\n" >> "$T/transcripts.tsv"
# T3: same gene+te as T1/T2 but a DIFFERENT chimera_type -- must stay its
# own row, proving the grouping key is 3-way (gene, te, chimera_type), not
# just (gene, te).
printf "T3\tMSTRG.3\tchr1\t+\t3\t1\t150\t60\t100\tTE_A\tL1PA2\tL1\tLINE\t2\tTE_A\tGENE1\t+\tGENE1\tyes\tte_exonized\n" >> "$T/transcripts.tsv"
# T4: matched_gene_id=="." (te_initiated_intergenic) -- not a real gene-TE
# candidate, must be excluded entirely.
printf "T4\tMSTRG.4\tchr1\t+\t2\t1\t100\t1\t50\tTE_B\tAluYa5\tAlu\tSINE\t1\tTE_B\t.\t.\t\tNA\tte_initiated_intergenic\n" >> "$T/transcripts.tsv"
# T5: independent singleton group -- sanity check unrelated groups don't merge.
printf "T5\tMSTRG.5\tchr2\t-\t2\t1\t100\t60\t100\tTE_C\tMER41\tERV1\tLTR\t2\tTE_C\tGENE2\t-\tGENE2\tyes\tte_terminated\n" >> "$T/transcripts.tsv"
# T6: te_id=="." -- classify_chimera_assembly.py itself never emits this
# (te_id is always resolved whenever chimera_type is set), but the
# aggregator's own filter must still catch it defensively.
printf "T6\tMSTRG.6\tchr2\t-\t2\t1\t100\t60\t100\t.\t.\t.\t.\t2\t.\tGENE2\t-\tGENE2\tyes\tte_terminated\n" >> "$T/transcripts.tsv"

printf "transcript_id\tS1\tS2\n" > "$T/counts.tsv"
printf "T1\t10\t20\n" >> "$T/counts.tsv"
printf "T2\t5\t8\n" >> "$T/counts.tsv"
printf "T3\t3\t4\n" >> "$T/counts.tsv"
printf "T4\t100\t200\n" >> "$T/counts.tsv"
printf "T5\t7\t9\n" >> "$T/counts.tsv"
printf "T6\t50\t60\n" >> "$T/counts.tsv"

if ! python3 workflow/scripts/aggregate_chimera_assembly_counts.py \
      --transcripts "$T/transcripts.tsv" --counts "$T/counts.tsv" \
      --out "$T/gene_te_chimera_counts.tsv" > "$T/aggregate.log" 2>&1; then
  echo "ERROR: aggregate_chimera_assembly_counts.py failed"; cat "$T/aggregate.log"; FAIL=1
else
  n_rows=$(($(wc -l < "$T/gene_te_chimera_counts.tsv") - 1))
  if [ "$n_rows" != "3" ]; then
    echo "ERROR: expected 3 output rows (T1+T2 summed, T3 its own row, T5 its own row -- T4/T6 excluded), got $n_rows"
    cat "$T/gene_te_chimera_counts.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE1:TE_A:te_initiated"{print $2"\t"$3}' "$T/gene_te_chimera_counts.tsv")
  if [ "$got" != "15	28" ]; then
    echo "ERROR: GENE1:TE_A:te_initiated should sum T1+T2 to 15/28, got '$got'"
    cat "$T/gene_te_chimera_counts.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE1:TE_A:te_exonized"{print $2"\t"$3}' "$T/gene_te_chimera_counts.tsv")
  if [ "$got" != "3	4" ]; then
    echo "ERROR: GENE1:TE_A:te_exonized (same gene+te as te_initiated, different chimera_type) must stay a separate row with T3's own counts 3/4, got '$got'"
    cat "$T/gene_te_chimera_counts.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE2:TE_C:te_terminated"{print $2"\t"$3}' "$T/gene_te_chimera_counts.tsv")
  if [ "$got" != "7	9" ]; then
    echo "ERROR: GENE2:TE_C:te_terminated singleton should be 7/9, got '$got'"
    cat "$T/gene_te_chimera_counts.tsv"; FAIL=1
  fi

  if grep -qE "100|200|	50	|	60$" "$T/gene_te_chimera_counts.tsv"; then
    echo "ERROR: T4 (matched_gene_id=.) and/or T6 (te_id=.) counts leaked into the output -- both should be excluded"
    cat "$T/gene_te_chimera_counts.tsv"; FAIL=1
  fi
fi

exit $FAIL
