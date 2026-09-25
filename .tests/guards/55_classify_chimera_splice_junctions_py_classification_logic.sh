#!/usr/bin/env bash
# Guard 55: classify_chimera_splice_junctions.py classification logic
#
# Run on its own:   .tests/guards/55_classify_chimera_splice_junctions_py_classification_logic.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- synthetic SJ.out.tab covering: a real te_to_gene/te_initiated call,
# --require-canonical filtering out a non-canonical junction, and
# --min-unique-reads filtering out a low-support one. The CI test fixture
# itself has ZERO real splice junctions (STAR: "Number of splices: Total |
# 0" -- the synthetic reads don't span any introns), so this is the only
# place this script's actual classification logic gets exercised.
printf "chr1\t1000\t1200\tGENE1\t.\t+\n" > "$T/genes.bed"
printf "chr1\t1000\t1200\tGENE1\t.\t+\n" > "$T/exons.bed"
printf "chr1\t500\t700\tTE_A\t.\t+\tL1\tLINE\tL1PA2\n" > "$T/te.bed"

# SJ.out.tab columns: chrom, intron_start, intron_end, strand(1=+), motif,
# annotated, unique_reads, multi_reads, overhang.
# J1: canonical, well-supported junction TE_A -> GENE1 (donor breakpoint at
#     699/700 overlaps TE_A [500,700); acceptor breakpoint at 1000/1001
#     overlaps GENE1's exon [1000,1200)) -- should be classified
#     te_to_gene / te_initiated (TE lies upstream of the gene on its own
#     strand).
# J2: same coordinates but motif 0 (non-canonical) -- must be dropped when
#     --require-canonical is passed, kept (and reported canonical=no) when
#     it is not.
# J3: canonical but only 1 unique read -- must be dropped by
#     --min-unique-reads 2.
printf "chr1\t700\t1000\t1\t1\t0\t5\t0\t20\n" > "$T/sj.tab"
printf "chr1\t700\t1000\t1\t0\t0\t5\t0\t20\n" >> "$T/sj.tab"
printf "chr1\t2700\t3000\t1\t1\t0\t1\t0\t20\n" >> "$T/sj.tab"

if ! python3 workflow/scripts/classify_chimera_splice_junctions.py \
      --sj "$T/sj.tab" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --sample S1 --breakpoint-tolerance 0 \
      --min-unique-reads 1 \
      --out "$T/lenient.tsv" > "$T/lenient.log" 2>&1; then
  echo "ERROR: classify_chimera_splice_junctions.py (lenient) failed"; cat "$T/lenient.log"; FAIL=1
else
  n_rows=$(($(wc -l < "$T/lenient.tsv") - 1))
  if [ "$n_rows" != "3" ]; then
    echo "ERROR: expected 3 rows with no filtering, got $n_rows"
    cat "$T/lenient.tsv"; FAIL=1
  fi
  j1_dir=$(awk -F'\t' 'NR==2{print $15}' "$T/lenient.tsv")
  j1_gene=$(awk -F'\t' 'NR==2{print $17}' "$T/lenient.tsv")
  j1_te=$(awk -F'\t' 'NR==2{print $19}' "$T/lenient.tsv")
  j1_type=$(awk -F'\t' 'NR==2{print $23}' "$T/lenient.tsv")
  j1_canon=$(awk -F'\t' 'NR==2{print $8}' "$T/lenient.tsv")
  if [ "$j1_dir" != "te_to_gene" ]; then
    echo "ERROR: J1 direction should be te_to_gene, got '$j1_dir'"; FAIL=1
  fi
  if [ "$j1_gene" != "GENE1" ] || [ "$j1_te" != "TE_A" ]; then
    echo "ERROR: J1 gene_id/te_id should be GENE1/TE_A, got '$j1_gene'/'$j1_te'"; FAIL=1
  fi
  if [ "$j1_type" != "te_initiated" ]; then
    echo "ERROR: J1 chimera_type should be te_initiated (TE upstream of GENE1 on + strand), got '$j1_type'"
    FAIL=1
  fi
  if [ "$j1_canon" != "yes" ]; then
    echo "ERROR: J1 canonical should be yes, got '$j1_canon'"; FAIL=1
  fi
  j2_canon=$(awk -F'\t' 'NR==3{print $8}' "$T/lenient.tsv")
  if [ "$j2_canon" != "no" ]; then
    echo "ERROR: J2 (motif 0) canonical should be no, got '$j2_canon'"; FAIL=1
  fi
fi

if ! python3 workflow/scripts/classify_chimera_splice_junctions.py \
      --sj "$T/sj.tab" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --sample S1 --breakpoint-tolerance 0 \
      --min-unique-reads 1 --require-canonical \
      --out "$T/strict.tsv" > "$T/strict.log" 2>&1; then
  echo "ERROR: classify_chimera_splice_junctions.py (--require-canonical) failed"; cat "$T/strict.log"; FAIL=1
else
  n_rows=$(($(wc -l < "$T/strict.tsv") - 1))
  if [ "$n_rows" != "2" ]; then
    echo "ERROR: --require-canonical should drop J2 (motif 0), expected 2 rows, got $n_rows"
    cat "$T/strict.tsv"; FAIL=1
  fi
fi

if ! python3 workflow/scripts/classify_chimera_splice_junctions.py \
      --sj "$T/sj.tab" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --sample S1 --breakpoint-tolerance 0 \
      --min-unique-reads 2 \
      --out "$T/min_reads.tsv" > "$T/min_reads.log" 2>&1; then
  echo "ERROR: classify_chimera_splice_junctions.py (--min-unique-reads 2) failed"; cat "$T/min_reads.log"; FAIL=1
else
  n_rows=$(($(wc -l < "$T/min_reads.tsv") - 1))
  if [ "$n_rows" != "2" ]; then
    echo "ERROR: --min-unique-reads 2 should drop J3 (1 read), expected 2 rows, got $n_rows"
    cat "$T/min_reads.tsv"; FAIL=1
  fi
fi

exit $FAIL
