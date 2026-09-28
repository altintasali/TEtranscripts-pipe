#!/usr/bin/env bash
# Guard 61: canonical_enrichment is a plain TSV data file, not its own
# report section
#
# BUG FIXED 2026: "Chimeric reads - splice-motif enrichment" used to be a
# 20-row Fisher table -- 16 of them per-sample repeats of the pooled
# result, not per-run QC (it characterises the classification method, not
# the cohort). The pooled result for each direction's primary comparison is
# now one sentence in the canonical-rate section's own description; the
# full per-sample+pooled table moved to
# results/chimera/chimeric_reads/canonical_enrichment.tsv.gz, a plain data
# file, not fed to MultiQC at all.
#
# Run on its own:   .tests/guards/61_canonical_enrichment_is_a_data_file_not_a_report_section.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

mkdir -p "$T/ce/qc"
# gene_to_te enriched over gene_to_gene (10/40=25% vs 3/60=5%); te_to_gene
# has no data at all, so that comparison must be silently skipped rather
# than crashing (0/0).
{ printf 'event_id\tsample\tdirection\tdirection_ambiguous\tcanonical\tchimera_type\tantisense_flag\tgene_strand_match\n'
  for i in $(seq 1 100); do
    if [ "$i" -le 10 ]; then printf 'e%s\ts1\tgene_to_te\tno\tyes\tte_initiated\tno\tyes\n' "$i"
    elif [ "$i" -le 40 ]; then printf 'e%s\ts1\tgene_to_te\tno\tno\tte_initiated\tno\tyes\n' "$i"
    elif [ "$i" -le 43 ]; then printf 'e%s\ts1\tgene_to_gene\tno\tyes\t.\t.\tNA\n' "$i"
    else printf 'e%s\ts1\tgene_to_gene\tno\tno\t.\t.\tNA\n' "$i"; fi
  done
} | gzip -c > "$T/ce/j.tsv.gz"
if ! python3 workflow/scripts/chimera_chimeric_reads_qc.py --table "$T/ce/j.tsv.gz" \
      --sample s1 --out "$T/ce/qc.tsv.gz" > "$T/ce/qc.log" 2>&1; then
  echo "ERROR: chimera_chimeric_reads_qc.py failed"; cat "$T/ce/qc.log"; FAIL=1
elif ! python3 workflow/scripts/chimera_chimeric_reads_qc_mqc.py \
      --tables "$T/ce/qc.tsv.gz" --samples s1 \
      --out "$T/ce/qc/chimera_chimeric_reads_qc_mqc.json" \
      --out-canonical "$T/ce/qc/canonical_rate_mqc.json" \
      --out-enrichment "$T/ce/canonical_enrichment.tsv.gz" \
      > "$T/ce/mqc.log" 2>&1; then
  echo "ERROR: chimera_chimeric_reads_qc_mqc.py failed"; cat "$T/ce/mqc.log"; FAIL=1
else
  # --- the TSV exists with the expected columns ------------------------
  if [ ! -f "$T/ce/canonical_enrichment.tsv.gz" ]; then
    echo "ERROR: canonical_enrichment.tsv.gz was not written"; FAIL=1
  else
    header=$(zcat "$T/ce/canonical_enrichment.tsv.gz" | head -1)
    want="comparison	sample	canonical	junctions	rate	comparator_rate	odds_ratio	p	q"
    if [ "$header" != "$want" ]; then
      echo "ERROR: canonical_enrichment.tsv.gz header wrong"
      echo "  want: $want"; echo "  got:  $header"; FAIL=1
    fi
    n_rows=$(($(zcat "$T/ce/canonical_enrichment.tsv.gz" | wc -l) - 1))
    if [ "$n_rows" -lt 1 ]; then
      echo "ERROR: canonical_enrichment.tsv.gz has no data rows"; FAIL=1
    fi
    if ! zcat "$T/ce/canonical_enrichment.tsv.gz" | grep -q "^gene_to_te vs gene_to_gene	pooled"; then
      echo "ERROR: expected pooled gene_to_te vs gene_to_gene row missing"
      zcat "$T/ce/canonical_enrichment.tsv.gz"; FAIL=1
    fi
    # the (a+b)==0 / (c+d)==0 skip must still apply -- no te_to_gene data
    # in this fixture, so that comparison contributes no rows at all.
    if zcat "$T/ce/canonical_enrichment.tsv.gz" | grep -q "te_to_gene vs te_to_te"; then
      echo "ERROR: an all-zero comparison (no te_to_gene data) leaked a row"
      zcat "$T/ce/canonical_enrichment.tsv.gz"; FAIL=1
    fi
  fi

  # --- the pooled result is folded into canonical-rate's description ---
  python3 - "$T/ce/qc/canonical_rate_mqc.json" <<'PY' || FAIL=1
import json, sys
d = json.load(open(sys.argv[1]))
desc = d.get("description", "")
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
check("Pooled across samples" in desc,
      f"canonical-rate description must carry the pooled summary sentence: {desc!r}")
check("canonical_enrichment.tsv.gz" in desc,
      f"canonical-rate description must point at the TSV: {desc!r}")
check("gene_to_te" in desc and "gene_to_gene" in desc,
      f"canonical-rate description must name the comparison it summarises: {desc!r}")
sys.exit(0 if ok else 1)
PY

  # --- no enrichment SECTION in the rendered report ---------------------
  if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
        -o "$T/ce/out" -n r "$T/ce/qc" > "$T/ce/render.log" 2>&1; then
    echo "ERROR: multiqc failed on the canonical-rate section"
    tail -30 "$T/ce/render.log"; FAIL=1
  else
    if grep -qi "splice-motif enrichment" "$T/ce/out/r.html"; then
      echo "ERROR: 'splice-motif enrichment' section text leaked back into the report"
      FAIL=1
    fi
    if grep -q "chimera_canonical_enrichment" "$T/ce/out/r.html"; then
      echo "ERROR: chimera_canonical_enrichment id leaked back into the report"
      FAIL=1
    fi
    # canonical_enrichment.tsv.gz sits outside qc/ (results/chimera/chimeric_reads/),
    # so it must not have been picked up as custom content even though it
    # was never copied into the scanned qc/ dir here.
    if grep -q "canonical_enrichment" "$T/ce/render.log"; then
      echo "ERROR: multiqc's own log mentions canonical_enrichment -- it should never scan for it"
      grep canonical_enrichment "$T/ce/render.log"; FAIL=1
    fi
  fi
fi

exit $FAIL
