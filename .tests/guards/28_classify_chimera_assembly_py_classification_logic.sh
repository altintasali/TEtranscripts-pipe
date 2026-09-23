#!/usr/bin/env bash
# Guard 28: classify_chimera_assembly.py classification logic
#
# Run on its own:   .tests/guards/28_classify_chimera_assembly_py_classification_logic.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- synthetic assembly covering all 6 output classes plus two
# negative controls: T5 (TE on the last exon, no earlier gene match --
# nothing to terminate, must be dropped, not just low-confidence) and
# T6 (TSS outside the TE that merely overlaps its first exon -- must be
# dropped under --require-tss-in-te even though it would be te_initiated
# without that flag).
printf "chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t20500\t20700\tGENE6\t.\t+\nchr1\t30500\t30700\tGENE7\t.\t+\n" > "$T/genes.bed"
printf "chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t20500\t20700\tGENE6\t.\t+\nchr1\t30500\t30700\tGENE7\t.\t+\n" > "$T/exons.bed"
printf "chr1\t500\t700\tTE_A\t.\t+\tL1\tLINE\tL1PA2\nchr1\t2600\t2800\tTE_B\t.\t+\tAluY\tSINE\tAluYa5\nchr1\t9000\t9100\tTE_C\t.\t+\tL2\tLINE\tL2a\nchr1\t1500\t1600\tTE_D\t.\t+\tERV1\tLTR\tMER41\nchr1\t8700\t8800\tTE_E\t.\t+\tL1\tLINE\tL1MA4\nchr1\t20150\t20250\tTE_F\t.\t+\tL1\tLINE\tL1MA5\nchr1\t30000\t30300\tTE_G\t.\t+\tERVL\tLTR\tMERVL\n" > "$T/te.bed"
# first_exons.bed: only T7's own first-exon interval is present, so T1..T6
# are unaffected by the annotated-promoter reclassification -- only T7
# (which would otherwise be a normal te_initiated call) should be
# downgraded to annotated_promoter_embedded_te.
printf "chr1\t30000\t30200\tSOME_TX\t.\t+\n" > "$T/first_exons.bed"
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
# T6: first exon's TE overlap does NOT cover the TSS (TE starts mid-exon) --
# a downstream exon matches GENE6, so this would be te_initiated under the
# loose (default) any-overlap behavior, but must be DROPPED entirely under
# --require-tss-in-te (falls through to te_terminated/te_exonized, neither
# of which apply here -- only 2 exons, so no internal exon to check either).
printf 'chr1\tStringTie\ttranscript\t20001\t20700\t.\t+\t.\ttranscript_id "T6"; gene_id "MSTRG.6";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t20001\t20200\t.\t+\t.\ttranscript_id "T6"; gene_id "MSTRG.6";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t20500\t20700\t.\t+\t.\ttranscript_id "T6"; gene_id "MSTRG.6";\n' >> "$T/stringtie.gtf"
# T7: TSS is inside TE_G (like T1), AND the first exon matches an annotated
# first exon (first_exons.bed) -- must be reclassified as
# annotated_promoter_embedded_te instead of te_initiated, even though a
# downstream exon also matches GENE7 (what would otherwise make it a normal
# te_initiated call).
printf 'chr1\tStringTie\ttranscript\t30001\t30700\t.\t+\t.\ttranscript_id "T7"; gene_id "MSTRG.7";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t30001\t30200\t.\t+\t.\ttranscript_id "T7"; gene_id "MSTRG.7";\n' >> "$T/stringtie.gtf"
printf 'chr1\tStringTie\texon\t30500\t30700\t.\t+\t.\ttranscript_id "T7"; gene_id "MSTRG.7";\n' >> "$T/stringtie.gtf"
if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --te "$T/te.bed" \
      --breakpoint-tolerance 0 --out "$T/candidates.tsv" > "$T/classify.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/classify.log"; FAIL=1
else
  declare -A want=( [T1]=te_initiated [T2]=te_terminated [T3]=te_exonized [T4]=unspliced_te_only
                     [T6]=te_initiated [T7]=annotated_promoter_embedded_te )
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
  t1_orient=$(awk -F'\t' '$1=="T1"{print $17}' "$T/candidates.tsv")
  if [ "$t1_orient" != "yes" ]; then
    echo "ERROR: T1's te_orientation_match should be 'yes' (both + strand), got '$t1_orient'"
    FAIL=1
  fi
fi

# --- --require-tss-in-te: T6's TE overlaps its first exon but not the TSS,
# so it must be dropped under this flag (unlike the default run above).
if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --te "$T/te.bed" --require-tss-in-te \
      --breakpoint-tolerance 0 --out "$T/candidates_strict.tsv" > "$T/classify_strict.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py --require-tss-in-te failed"
  cat "$T/classify_strict.log"; FAIL=1
elif grep -qP "^T6\t" "$T/candidates_strict.tsv"; then
  echo "ERROR: T6 (TE overlaps first exon but not the TSS) should be dropped under --require-tss-in-te"
  cat "$T/candidates_strict.tsv"; FAIL=1
fi

exit $FAIL
