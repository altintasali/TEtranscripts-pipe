#!/usr/bin/env bash
# Guard 58: chimera_splice_junctions_counts.py's counts_matrix/cpm_matrix
# (and te-gene-junctions.tsv) are restricted to gene<->TE events (direction
# gene_to_te / te_to_gene); gene_to_gene/other events stay in the full
# all_events catalog but must NOT leak into the matrix or the QC view it
# feeds. Also pins the CPM denominator: it must be computed from EVERY
# event's reads (the sample's whole splice-junction library), not just the
# TE-restricted rows the matrix reports -- otherwise CPM silently stops
# meaning "per million reads in this sample" and instead means "per million
# TE-junction reads", which is a different, incompatible number.
#
# Measured 2026 (see chimera_splice_junctions_counts.py's module docstring):
# on real data, gene<->TE events were ~4% of all classified junctions (76,430
# of ~2M across 12 real mouse samples); on the unrestricted matrix, DESeq2's
# vst() failed outright on a synthetic worst-case cohort ("every gene
# contains at least one zero") and silently fell back to log2 -- restricting
# the matrix to gene<->TE events fixed both the QC-quality regression and
# the ~13x oversized matrix, without needing to touch all_events itself
# (measured cheap even unrestricted: 963 MB / 3.5 GB peak RSS at 12 / 80
# samples respectively, well under the 8 GB resources.yaml budget).
#
# Run on its own:   .tests/guards/58_chimera_splice_junctions_counts_py_matrix_is_te_restricted.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

header="event_id\tsample\tchrom\tintron_start\tintron_end\tstrand\tmotif\tcanonical\tannotated\tunique_reads\tmulti_reads\toverhang\tdonor_hits\tacceptor_hits\tdirection\tdirection_ambiguous\tgene_id\tgene_strand\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\tte_initiated_detail\tantisense_flag\tlibrary_strand\ttranscript_strand\tgene_strand_match"

# S1: EVT1 (gene_to_te, 10 reads), EVT2 (te_to_gene, 20 reads), EVT3
# (gene_to_gene, 900 reads -- the overwhelming majority of real junctions),
# EVT4 (other, 70 reads). Total = 1000, chosen so CPM math is exact.
printf "$header\n" > "$T/S1.tsv"
printf "EVT1\tS1\tchr1\t100\t200\t+\t1\tyes\t1\t10\t0\t20\t.\t.\tgene_to_te\tno\tGENE1\t+\tTE1\t.\t.\t.\tte_terminated\t.\t.\tNA\tNA\tNA\n" >> "$T/S1.tsv"
printf "EVT2\tS1\tchr1\t300\t400\t+\t1\tyes\t1\t20\t0\t20\t.\t.\tte_to_gene\tno\tGENE1\t+\tTE2\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\n" >> "$T/S1.tsv"
printf "EVT3\tS1\tchr1\t500\t600\t+\t1\tyes\t1\t900\t0\t20\t.\t.\tgene_to_gene\tno\tGENE1\t+\t.\t.\t.\t.\t.\t.\t.\tNA\tNA\tNA\n" >> "$T/S1.tsv"
printf "EVT4\tS1\tchr1\t700\t800\t+\t1\tyes\t1\t70\t0\t20\t.\t.\tother\tno\t.\t.\t.\t.\t.\t.\t.\t.\t.\tNA\tNA\tNA\n" >> "$T/S1.tsv"

# S2: EVT1 (5 reads), EVT2 (15 reads), EVT3 (770 reads), EVT4 (210 reads).
# Total = 1000 again.
printf "$header\n" > "$T/S2.tsv"
printf "EVT1\tS2\tchr1\t100\t200\t+\t1\tyes\t1\t5\t0\t20\t.\t.\tgene_to_te\tno\tGENE1\t+\tTE1\t.\t.\t.\tte_terminated\t.\t.\tNA\tNA\tNA\n" >> "$T/S2.tsv"
printf "EVT2\tS2\tchr1\t300\t400\t+\t1\tyes\t1\t15\t0\t20\t.\t.\tte_to_gene\tno\tGENE1\t+\tTE2\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\n" >> "$T/S2.tsv"
printf "EVT3\tS2\tchr1\t500\t600\t+\t1\tyes\t1\t770\t0\t20\t.\t.\tgene_to_gene\tno\tGENE1\t+\t.\t.\t.\t.\t.\t.\t.\tNA\tNA\tNA\n" >> "$T/S2.tsv"
printf "EVT4\tS2\tchr1\t700\t800\t+\t1\tyes\t1\t210\t0\t20\t.\t.\tother\tno\t.\t.\t.\t.\t.\t.\t.\t.\t.\tNA\tNA\tNA\n" >> "$T/S2.tsv"

if ! python3 workflow/scripts/chimera_splice_junctions_counts.py \
      --tables "$T/S1.tsv" "$T/S2.tsv" --sample-names S1 S2 \
      --out-events "$T/all_events.tsv" \
      --out-counts "$T/counts_matrix.tsv" \
      --out-cpm "$T/cpm_matrix.tsv" \
      --out-te-events "$T/te_gene_junctions.tsv" \
      > "$T/merge.log" 2>&1; then
  echo "ERROR: chimera_splice_junctions_counts.py failed"; cat "$T/merge.log"; FAIL=1
else
  # all_events: every direction class present.
  n_all=$(($(wc -l < "$T/all_events.tsv") - 1))
  if [ "$n_all" != "4" ]; then
    echo "ERROR: all_events.tsv should keep all 4 event classes, got $n_all rows"
    cat "$T/all_events.tsv"; FAIL=1
  fi

  # counts_matrix / cpm_matrix / te-gene-junctions: ONLY EVT1/EVT2 (gene<->TE).
  n_counts=$(($(wc -l < "$T/counts_matrix.tsv") - 1))
  if [ "$n_counts" != "2" ]; then
    echo "ERROR: counts_matrix.tsv should have exactly 2 rows (EVT1, EVT2 -- gene<->TE only), got $n_counts"
    cat "$T/counts_matrix.tsv"; FAIL=1
  fi
  if grep -qE "^EVT3|^EVT4" "$T/counts_matrix.tsv" "$T/cpm_matrix.tsv" "$T/te_gene_junctions.tsv"; then
    echo "ERROR: gene_to_gene/other events (EVT3/EVT4) leaked into the TE-restricted matrix or te-gene-junctions catalog"
    FAIL=1
  fi
  if ! grep -q "^EVT1" "$T/counts_matrix.tsv" || ! grep -q "^EVT2" "$T/counts_matrix.tsv"; then
    echo "ERROR: EVT1/EVT2 (gene<->TE events) missing from counts_matrix.tsv"
    cat "$T/counts_matrix.tsv"; FAIL=1
  fi

  # counts_matrix values themselves (not just row presence).
  got=$(awk -F'\t' '$1=="EVT1"{print $2"\t"$3}' "$T/counts_matrix.tsv")
  if [ "$got" != "10	5" ]; then
    echo "ERROR: EVT1 counts should be S1=10, S2=5, got '$got'"; FAIL=1
  fi

  # CPM denominator must be the FULL per-sample total (1000 for both S1 and
  # S2 here, by construction), not the TE-subset total (30 for S1, 20 for
  # S2) -- that's the whole point of this guard. EVT1 in S1: 10/1000*1e6 =
  # 10000.000. If the denominator were wrongly TE-subset-only (10+20=30),
  # this would instead read 333333.333.
  got=$(awk -F'\t' '$1=="EVT1"{print $2}' "$T/cpm_matrix.tsv")
  if [ "$got" != "10000.000" ]; then
    echo "ERROR: EVT1 S1 CPM should be 10000.000 (denominator = whole-library total 1000, not TE-subset total 30), got '$got' -- CPM denominator is no longer using all events"
    cat "$T/cpm_matrix.tsv"; FAIL=1
  fi
  got=$(awk -F'\t' '$1=="EVT2"{print $3}' "$T/cpm_matrix.tsv")
  if [ "$got" != "15000.000" ]; then
    echo "ERROR: EVT2 S2 CPM should be 15000.000 (denominator = whole-library total 1000, not TE-subset total 20), got '$got'"
    cat "$T/cpm_matrix.tsv"; FAIL=1
  fi
fi

exit $FAIL
