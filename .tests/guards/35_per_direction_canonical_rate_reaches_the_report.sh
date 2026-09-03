#!/usr/bin/env bash
# Guard 35: per-direction canonical rate reaches the report
#
# Run on its own:   .tests/guards/35_per_direction_canonical_rate_reaches_the_report.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# "canonical" = STAR reported any recognised splice motif. Real
# introns are ~100% canonical, so this is the main signal-vs-
# artifact discriminator -- and it was previously reported only
# GLOBALLY, which hid that the gene-TE classes sit several-fold
# above the unannotated background. Fixture: gene_to_te 10/40
# canonical (25%), other 4/60 (6.7%).
mkdir -p "$T/canon/qc"
{ printf 'event_id\tsample\tdirection\tdirection_ambiguous\tcanonical\tchimera_type\tantisense_flag\tgene_strand_match\n'
  for i in $(seq 1 100); do
    if [ "$i" -le 10 ]; then printf 'e%s\ts\tgene_to_te\tno\tyes\tte_initiated\tno\tyes\n' "$i"
    elif [ "$i" -le 40 ]; then printf 'e%s\ts\tgene_to_te\tno\tno\tte_initiated\tno\tyes\n' "$i"
    elif [ "$i" -le 44 ]; then printf 'e%s\ts\tother\tno\tyes\t.\t.\tNA\n' "$i"
    else printf 'e%s\ts\tother\tno\tno\t.\t.\tNA\n' "$i"; fi
  done
} | gzip -c > "$T/canon/j.tsv.gz"
if ! python3 workflow/scripts/chimera_reads_qc.py --table "$T/canon/j.tsv.gz" \
      --sample s --out "$T/canon/qc.tsv.gz" > "$T/canon/qc.log" 2>&1; then
  echo "ERROR: chimera_reads_qc.py failed"; cat "$T/canon/qc.log"; FAIL=1
elif ! zcat "$T/canon/qc.tsv.gz" | grep -qx "canonical_gene_to_te	10"; then
  echo "ERROR: canonical_gene_to_te should be 10"
  zcat "$T/canon/qc.tsv.gz" | grep canonical; FAIL=1
elif ! python3 workflow/scripts/chimera_reads_qc_mqc.py --tables "$T/canon/qc.tsv.gz" \
      --samples s --out "$T/canon/qc/chimera_reads_qc_mqc.json" \
      --out-canonical "$T/canon/qc/canonical_rate_mqc.json" > "$T/canon/mqc.log" 2>&1; then
  echo "ERROR: chimera_reads_qc_mqc.py failed"; cat "$T/canon/mqc.log"; FAIL=1
else
  # Dataset ORDER is load-bearing: raw counts first, the derived rate second.
  # It used to be the other way round, which put the derived number in front
  # of the measurement and made the toggle read backwards.
  # 10/40 = 25.0% for gene_to_te, 4/60 = 6.7% for other.
  if ! python3 - "$T/canon/qc/canonical_rate_mqc.json" <<'PY2'
import json, sys
d = json.load(open(sys.argv[1]))
ok = True
labels = [x["name"] for x in d["pconfig"]["data_labels"]]
if labels != ["Canonical junctions", "% canonical"]:
    print(f"ERROR: datasets must be counts then rate; got {labels}"); ok = False
counts, rate = d["data"][0]["s"], d["data"][1]["s"]
if counts.get("gene_to_te") != 10:
    print(f"ERROR: first dataset must be raw counts; got {counts}"); ok = False
if not (abs(rate["gene_to_te"] - 25.0) < 0.05 and abs(rate["other"] - 6.7) < 0.05):
    print(f"ERROR: canonical percentages wrong: {rate}"); ok = False
# a rate is not a share of a total, so cpswitch must stay off -- its
# percentage would be a different number
if d["pconfig"].get("cpswitch") is not False:
    print("ERROR: cpswitch must be off for a rate plot"); ok = False
# and the 0-100 cap belongs to the percentage view only, or it clips counts
if d["pconfig"].get("ymax") is not None:
    print("ERROR: ymax at plot level clips the counts view"); ok = False
if d["pconfig"]["data_labels"][1].get("ymax") != 100:
    print("ERROR: the percentage dataset must cap at 100"); ok = False
sys.exit(0 if ok else 1)
PY2
  then
    cat "$T/canon/qc/canonical_rate_mqc.json"; FAIL=1
  fi
  if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
        -o "$T/canon/mqc" -n r "$T/canon/qc" > "$T/canon/render.log" 2>&1; then
    echo "ERROR: multiqc failed on the canonical section"; tail -30 "$T/canon/render.log"; FAIL=1
  elif ! grep -q "custom_content | chimera_canonical_rate: Found" "$T/canon/render.log"; then
    echo "ERROR: canonical-rate section not rendered"
    grep custom_content "$T/canon/render.log" | head; FAIL=1
  fi
fi

# ...and the ALL-ZERO case, which is what the synthetic fixture
# actually produces (its chimeric reads are random sequence, so no
# junction carries a splice motif). An all-zero bar document makes
# MultiQC raise "No datasets to plot" and fails the whole report,
# so it must degrade to html. The fixture above has canonical hits
# and therefore could not catch this.
mkdir -p "$T/canon0/qc"
{ printf 'event_id\tsample\tdirection\tdirection_ambiguous\tcanonical\tchimera_type\tantisense_flag\tgene_strand_match\n'
  for i in $(seq 1 30); do printf 'e%s\ts\tgene_to_te\tno\tno\tte_initiated\tno\tyes\n' "$i"; done
} | gzip -c > "$T/canon0/j.tsv.gz"
python3 workflow/scripts/chimera_reads_qc.py --table "$T/canon0/j.tsv.gz" \
  --sample s --out "$T/canon0/qc.tsv.gz" > /dev/null 2>&1
if ! python3 workflow/scripts/chimera_reads_qc_mqc.py --tables "$T/canon0/qc.tsv.gz" \
      --samples s --out "$T/canon0/qc/chimera_reads_qc_mqc.json" \
      --out-canonical "$T/canon0/qc/canonical_rate_mqc.json" > "$T/canon0/mqc.log" 2>&1; then
  echo "ERROR: chimera_reads_qc_mqc.py failed on all-zero canonical"; cat "$T/canon0/mqc.log"; FAIL=1
elif ! python3 -c "import json,sys; sys.exit(0 if json.load(open('$T/canon0/qc/canonical_rate_mqc.json'))['plot_type']=='html' else 1)"; then
  echo "ERROR: all-zero canonical must degrade to plot_type html"; FAIL=1
elif ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/canon0/mqc_out" -n r "$T/canon0/qc" > "$T/canon0/render.log" 2>&1; then
  echo "ERROR: multiqc failed on an all-zero canonical rate"
  tail -30 "$T/canon0/render.log"; FAIL=1
fi

exit $FAIL
