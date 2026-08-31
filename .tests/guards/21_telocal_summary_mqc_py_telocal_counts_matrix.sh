#!/usr/bin/env bash
# shellcheck disable=SC2010  # ls|grep over known-safe generated filenames
# Guard 21: telocal_summary_mqc.py + telocal counts matrix
#
# Run on its own:   .tests/guards/21_telocal_summary_mqc_py_telocal_counts_matrix.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- telocal_summary_mqc.py must render in MultiQC ------------------
printf 'gene/TE\t/path/a.bam\nENSG0001\t80\nchr1:100:200(L1PA2:+):L1PA2:L1:LINE\t5\n' > "$T/tl_a.cntTable"
printf 'gene/TE\t/path/b.bam\nENSG0001\t0\nchr1:100:200(L1PA2:+):L1PA2:L1:LINE\t0\n' > "$T/tl_b.cntTable"
mkdir -p "$T/custom/telocal/qc"
if ! python workflow/scripts/telocal_summary_mqc.py \
      --tables "$T/tl_a.cntTable" "$T/tl_b.cntTable" --samples a b \
      --out-assignment "$T/custom/telocal/qc/telocal_assignment_mqc.json" \
      --out-class "$T/custom/telocal/qc/telocal_te_class_mqc.json" > "$T/tlsum.log" 2>&1; then
  echo "ERROR: telocal_summary_mqc.py failed"; cat "$T/tlsum.log"; FAIL=1
fi

# --- tecount_counts.py --key-style telocal classifies by shape ------
# TE-locus keys (>= 3 colon fields) are kept with --feature-class TE;
# gene keys (no colons) are dropped.
if ! python workflow/scripts/tecount_counts.py \
      --tables "$T/tl_a.cntTable" "$T/tl_b.cntTable" --sample-names a b \
      --key-style telocal --feature-class TE \
      --out-counts "$T/telocal_counts_matrix.tsv" > "$T/tlcounts.log" 2>&1; then
  echo "ERROR: tecount_counts.py --key-style telocal failed"; cat "$T/tlcounts.log"; FAIL=1
elif ! grep -q "^chr1:100:200(L1PA2:+):L1PA2:L1:LINE" "$T/telocal_counts_matrix.tsv" \
   || grep -q "^ENSG0001" "$T/telocal_counts_matrix.tsv"; then
  echo "ERROR: telocal counts matrix kept wrong features (TE-only filter)"
  cat "$T/telocal_counts_matrix.tsv"; FAIL=1
fi

if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/mqc_tl" -n report "$T/custom" > "$T/mqc_tl.log" 2>&1; then
  echo "ERROR: multiqc run with telocal failed"; tail -40 "$T/mqc_tl.log"; FAIL=1
else
  for cid in telocal_assignment telocal_te_class; do
    if ! grep -q "custom_content | $cid: Found" "$T/mqc_tl.log" \
       || [ ! -f "$T/mqc_tl/report_data/multiqc_${cid}_plot.txt" ]; then
      echo "ERROR: section '$cid' not rendered by MultiQC"
      grep "custom_content" "$T/mqc_tl.log" | head -20
      ls "$T/mqc_tl/report_data" | grep "${cid}" || true
      FAIL=1
    fi
  done
fi

exit $FAIL
