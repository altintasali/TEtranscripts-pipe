#!/usr/bin/env bash
# Guard 54: aggregate_chimera_assembly_counts.py groups by
# (matched_gene_id, te_id, chimera_type), not by transcript_id or by
# (gene, te) alone, drops non-candidate rows, and the companion annotation
# table carries the right gene symbol/locus/TE annotation/traceability.
#
# Run on its own:   .tests/guards/54_aggregate_chimera_assembly_counts_py_groups_by_gene_te_chimera_type.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- transcripts.tsv.gz-shaped fixture (classify_chimera_assembly.py's own
# header). Only transcript_id(1), te_id(10), te_subfamily/family/class
# (11-13), matched_gene_id(16), and chimera_type(20) matter here; the rest
# are cheap placeholders.
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
# T5: independent singleton group -- sanity check unrelated groups don't
# merge, and its gene has no gene_names entry (must fall back to gene_id).
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

# GENE2 deliberately has no row here, to pin the gene_symbol fallback to
# gene_id when a gene isn't in --gene-names.
printf "gene_id\tgene_name\n" > "$T/gene_names.tsv"
printf "GENE1\tSymbolOne\n" >> "$T/gene_names.tsv"

# genes.bed / te.bed: chrom start end id score strand [family class
# subfamily] -- 0-based start, no header (annotation_to_bed.py's
# convention). GENE2/TE_C deliberately absent, to pin the "." fallback when
# a locus isn't in the BED.
printf "chr1\t999\t2000\tGENE1\t.\t+\n" > "$T/genes.bed"
printf "chr1\t2099\t2400\tTE_A\t.\t+\tL1\tLINE\tL1PA2\n" > "$T/te.bed"

if ! python3 workflow/scripts/aggregate_chimera_assembly_counts.py \
      --transcripts "$T/transcripts.tsv" --counts "$T/counts.tsv" \
      --gene-names "$T/gene_names.tsv" \
      --genes-bed "$T/genes.bed" --te-bed "$T/te.bed" \
      --out-counts "$T/counts_out.tsv" --out-annotation "$T/annotation_out.tsv" \
      > "$T/aggregate.log" 2>&1; then
  echo "ERROR: aggregate_chimera_assembly_counts.py failed"; cat "$T/aggregate.log"; FAIL=1
else
  n_rows=$(($(wc -l < "$T/counts_out.tsv") - 1))
  if [ "$n_rows" != "3" ]; then
    echo "ERROR: expected 3 counts-matrix rows (T1+T2 summed, T3 its own row, T5 its own row -- T4/T6 excluded), got $n_rows"
    cat "$T/counts_out.tsv"; FAIL=1
  fi
  n_annot_rows=$(($(wc -l < "$T/annotation_out.tsv") - 1))
  if [ "$n_annot_rows" != "3" ]; then
    echo "ERROR: expected 3 annotation rows, got $n_annot_rows"
    cat "$T/annotation_out.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE1:TE_A:te_initiated"{print $2"\t"$3}' "$T/counts_out.tsv")
  if [ "$got" != "15	28" ]; then
    echo "ERROR: GENE1:TE_A:te_initiated should sum T1+T2 to 15/28, got '$got'"
    cat "$T/counts_out.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE1:TE_A:te_exonized"{print $2"\t"$3}' "$T/counts_out.tsv")
  if [ "$got" != "3	4" ]; then
    echo "ERROR: GENE1:TE_A:te_exonized (same gene+te as te_initiated, different chimera_type) must stay a separate row with T3's own counts 3/4, got '$got'"
    cat "$T/counts_out.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE2:TE_C:te_terminated"{print $2"\t"$3}' "$T/counts_out.tsv")
  if [ "$got" != "7	9" ]; then
    echo "ERROR: GENE2:TE_C:te_terminated singleton should be 7/9, got '$got'"
    cat "$T/counts_out.tsv"; FAIL=1
  fi

  if grep -qE "100|200|	50	|	60$" "$T/counts_out.tsv"; then
    echo "ERROR: T4 (matched_gene_id=.) and/or T6 (te_id=.) counts leaked into the counts matrix -- both should be excluded"
    cat "$T/counts_out.tsv"; FAIL=1
  fi

  # --- annotation table: gene_symbol resolved (GENE1 -> SymbolOne),
  # gene_symbol fallback (GENE2, no gene_names entry -> GENE2 itself),
  # locus joins (present for GENE1/TE_A, "." fallback for GENE2/TE_C),
  # te_subfamily/family/class from the TE GTF (not StringTie), n_transcripts
  # (2 for the summed te_initiated group), and traceability back to the
  # source transcript_ids.
  got=$(awk -F'\t' '$1=="GENE1:TE_A:te_initiated"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9\
        "\t"$10"\t"$11"\t"$12}' "$T/annotation_out.tsv")
  want="GENE1	SymbolOne	chr1:1000-2000	TE_A	chr1:2100-2400	L1PA2	L1	LINE	te_initiated	2	T1,T2"
  if [ "$got" != "$want" ]; then
    echo "ERROR: GENE1:TE_A:te_initiated annotation row wrong."
    echo "  want: $want"
    echo "  got:  $got"
    cat "$T/annotation_out.tsv"; FAIL=1
  fi

  got=$(awk -F'\t' '$1=="GENE2:TE_C:te_terminated"{print $2"\t"$3"\t"$4"\t"$5}' "$T/annotation_out.tsv")
  want="GENE2	GENE2	.	TE_C"
  if [ "$got" != "$want" ]; then
    echo "ERROR: GENE2:TE_C:te_terminated annotation row should fall back to gene_symbol=GENE2 (no gene_names entry) and gene_locus=. (not in genes.bed), want '$want' got '$got'"
    cat "$T/annotation_out.tsv"; FAIL=1
  fi

  if grep -qE "(^|[,	])T4([,	]|$)|(^|[,	])T6([,	]|$)" "$T/annotation_out.tsv"; then
    echo "ERROR: T4/T6 (excluded transcripts) leaked into the annotation table's assembly_transcript_ids"
    cat "$T/annotation_out.tsv"; FAIL=1
  fi
fi

exit $FAIL
