#!/usr/bin/env bash
# Guard 64: classify_chimera_chimeric_reads.py classification logic
#
# Mirrors guard 55 (the splice-junctions screen's own classification guard):
# chimera_type here used to be decided from where the TE's genomic SPAN sits
# relative to the gene's overall span, which ignored the junction's own
# DIRECTION. Fixed to mirror classify_chimera_splice_junctions.py exactly,
# via the shared chimera_exon_context.py helpers -- see either script's
# module docstring for the full reasoning.
#
# Run on its own:   .tests/guards/64_classify_chimera_chimeric_reads_py_classification_logic.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# genes.bed / exons.bed: 6 cols (chrom, start, end, gene_id, score, strand).
# te.bed: 9 cols (chrom, start, end, te_id, score, strand, family, class,
# subfamily) -- same layout as guard 55's TE_A row.
#
# Seven cases on non-overlapping 20,000 bp blocks of chr1 (plus one TE on
# chr2 for the trans case), one gene<->TE junction each:
#   A  (+) te_to_gene,  TE upstream of gene's only exon      -> te_initiated / upstream
#   B  (+) te_to_gene,  intronic TE into a 2nd exon,          -> te_initiated / internal
#          skipping the gene's 1st exon (REGRESSION: old code    (old code: te_exonized)
#          called this te_exonized because the TE's span sits
#          inside the gene's overall genomic span)
#   B2 (-) same as B, mirrored to the minus strand             -> te_initiated / internal
#   C  (+) gene_to_te,  no exon downstream of the donor exon  -> te_terminated
#   D  (+) gene_to_te,  an exon downstream of the donor exon  -> te_exonized
#   E  (+) gene_to_te,  TE genomically UPSTREAM of the gene   -> te_terminated, NEVER
#          (REGRESSION: old code called this te_initiated         te_initiated
#          purely from TE-vs-gene span position, ignoring
#          that the junction's direction is gene_to_te)
#   F  (+) gene_to_te,  donor/acceptor on different chroms    -> "." (trans)
mkdir -p "$T/cr"
{
  printf 'chr1\t1999\t2300\tGENE_A\t.\t+\n'
  printf 'chr1\t21999\t25300\tGENE_B\t.\t+\n'
  printf 'chr1\t41999\t45300\tGENE_B2\t.\t-\n'
  printf 'chr1\t61999\t62300\tGENE_C\t.\t+\n'
  printf 'chr1\t81999\t85300\tGENE_D\t.\t+\n'
  printf 'chr1\t105000\t106000\tGENE_E\t.\t+\n'
  printf 'chr1\t121999\t122300\tGENE_F\t.\t+\n'
} > "$T/cr/genes.bed"

{
  printf 'chr1\t1999\t2300\tGENE_A\t.\t+\n'
  printf 'chr1\t21999\t22300\tGENE_B\t.\t+\n'
  printf 'chr1\t24999\t25300\tGENE_B\t.\t+\n'
  printf 'chr1\t44999\t45300\tGENE_B2\t.\t-\n'
  printf 'chr1\t41999\t42300\tGENE_B2\t.\t-\n'
  printf 'chr1\t61999\t62300\tGENE_C\t.\t+\n'
  printf 'chr1\t81999\t82300\tGENE_D\t.\t+\n'
  printf 'chr1\t84999\t85300\tGENE_D\t.\t+\n'
  printf 'chr1\t105000\t106000\tGENE_E\t.\t+\n'
  printf 'chr1\t121999\t122300\tGENE_F\t.\t+\n'
} > "$T/cr/exons.bed"

{
  printf 'chr1\t100\t500\tTE_A\t.\t+\tL1\tLINE\tL1PA2\n'
  printf 'chr1\t22999\t23500\tTE_B\t.\t+\tL1\tLINE\tL1PA2\n'
  printf 'chr1\t42999\t43500\tTE_B2\t.\t-\tL1\tLINE\tL1PA2\n'
  printf 'chr1\t63000\t63500\tTE_C\t.\t+\tL1\tLINE\tL1PA2\n'
  printf 'chr1\t83000\t83500\tTE_D\t.\t+\tL1\tLINE\tL1PA2\n'
  printf 'chr1\t100100\t100500\tTE_E\t.\t+\tL1\tLINE\tL1PA2\n'
  printf 'chr2\t100\t500\tTE_F\t.\t+\tL1\tLINE\tL1PA2\n'
} > "$T/cr/te.bed"

# Chimeric.out.junction columns: donor_chrom, donor_bp, donor_strand,
# acceptor_chrom, acceptor_bp, acceptor_strand, junction_type, repeat_flag.
# Column 2 is the first base of the DONOR'S INTRON and column 5 the last
# base of the ACCEPTOR'S INTRON (see donor_locus/acceptor_locus and guard
# 34) -- breakpoints below are placed exactly on the relevant feature edges.
{
  printf 'chr1\t501\t+\tchr1\t1999\t+\t1\t0\n'       # A
  printf 'chr1\t23501\t+\tchr1\t24999\t+\t1\t0\n'    # B
  printf 'chr1\t43499\t-\tchr1\t42301\t-\t1\t0\n'    # B2
  printf 'chr1\t62301\t+\tchr1\t63000\t+\t1\t0\n'    # C
  printf 'chr1\t82301\t+\tchr1\t83000\t+\t1\t0\n'    # D
  printf 'chr1\t106001\t+\tchr1\t100100\t+\t1\t0\n'  # E
  printf 'chr1\t122301\t+\tchr2\t100\t+\t1\t0\n'     # F
} > "$T/cr/j.junction"

if ! python3 workflow/scripts/classify_chimera_chimeric_reads.py \
      --junctions "$T/cr/j.junction" --genes "$T/cr/genes.bed" \
      --exons "$T/cr/exons.bed" --te "$T/cr/te.bed" --sample S1 \
      --breakpoint-tolerance 0 --library-strandedness no \
      --out "$T/cr/out.tsv" > "$T/cr/parse.log" 2>&1; then
  echo "ERROR: classify_chimera_chimeric_reads.py failed"; cat "$T/cr/parse.log"; FAIL=1
else
  # field(gene_id, column_name) -> value of that column on the one row
  # whose gene_id matches, via a header-name lookup (not a hardcoded
  # position) -- same idiom as guard 34.
  field() {
    awk -F'\t' -v gid="$1" -v col="$2" \
      'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["gene_id"]==gid{print $h[col]}' \
      "$T/cr/out.tsv"
  }

  check() {
    # check <label> <gene_id> <column> <expected>
    got=$(field "$2" "$3")
    if [ "$got" != "$4" ]; then
      echo "ERROR: $1 -- expected $3=$4, got '$got'"
      FAIL=1
    fi
  }

  # A: te_to_gene, TE upstream of the gene's only exon -> te_initiated / upstream
  check "A direction" GENE_A direction te_to_gene
  check "A chimera_type" GENE_A chimera_type te_initiated
  check "A te_initiated_detail" GENE_A te_initiated_detail upstream

  # B: intronic TE splicing past GENE_B's 1st exon into its 2nd -> te_initiated /
  # internal. Old span-based code called this te_exonized (TE's coordinates
  # fall inside the gene's overall span) -- assert the regression explicitly.
  check "B direction" GENE_B direction te_to_gene
  check "B chimera_type (regression: must not be te_exonized)" GENE_B chimera_type te_initiated
  check "B te_initiated_detail (regression: must not be upstream)" GENE_B te_initiated_detail internal

  # B2: same intronic-promoter case, minus strand.
  check "B2 direction" GENE_B2 direction te_to_gene
  check "B2 chimera_type" GENE_B2 chimera_type te_initiated
  check "B2 te_initiated_detail" GENE_B2 te_initiated_detail internal

  # C: gene_to_te, nothing annotated downstream of the donor exon -> te_terminated
  check "C direction" GENE_C direction gene_to_te
  check "C chimera_type" GENE_C chimera_type te_terminated
  check "C te_initiated_detail (gene_to_te always '.')" GENE_C te_initiated_detail .

  # D: gene_to_te, an exon annotated downstream of the donor exon -> te_exonized
  check "D direction" GENE_D direction gene_to_te
  check "D chimera_type" GENE_D chimera_type te_exonized

  # E: gene_to_te into a TE that sits genomically UPSTREAM of the gene. Old
  # span-based code called this te_initiated (TE span before gene span,
  # direction ignored entirely) -- assert the regression explicitly: this
  # must type like any other gene_to_te event (no downstream exon here), and
  # must never be te_initiated.
  check "E direction" GENE_E direction gene_to_te
  check "E chimera_type (regression: must not be te_initiated)" GENE_E chimera_type te_terminated

  # F: donor and acceptor on different chromosomes -> chimera_type stays "."
  check "F direction" GENE_F direction gene_to_te
  check "F chimera_type (trans event)" GENE_F chimera_type .
fi

exit $FAIL
