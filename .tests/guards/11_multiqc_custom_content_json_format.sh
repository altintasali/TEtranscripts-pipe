#!/usr/bin/env bash
# shellcheck disable=SC2010  # ls|grep over known-safe generated filenames
# Guard 11: MultiQC custom-content JSON format
#
# Run on its own:   .tests/guards/11_multiqc_custom_content_json_format.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- the barplot emitters must render in the installed MultiQC ------
# MultiQC's custom_content module only reads a top-level "data" key;
# the legacy "samples"/"datasets" shape yields "No data found" and the
# section is silently dropped (regressed in v0.6.1). Run both scripts
# on synthetic fixtures and confirm MultiQC renders all four sections.
printf 'metric\tvalue\nsample\ta\nevents_total\t10\ndirection_gene_to_te\t2\ndirection_te_to_gene\t1\ndirection_gene_to_gene\t4\ndirection_te_to_te\t0\ndirection_gene_to_other\t1\ndirection_other_to_gene\t2\ndirection_te_to_other\t0\ndirection_other_to_te\t1\ndirection_other\t0\n' > "$T/jqc_a.tsv"
printf 'metric\tvalue\nsample\tb\nevents_total\t0\ndirection_gene_to_te\t0\ndirection_te_to_gene\t0\ndirection_gene_to_gene\t0\ndirection_te_to_te\t0\ndirection_gene_to_other\t0\ndirection_other_to_gene\t0\ndirection_te_to_other\t0\ndirection_other_to_te\t0\ndirection_other\t0\n' > "$T/jqc_b.tsv"
printf 'gene/TE\t/path/a.bam\nENSG0001\t80\nL1PA_0001:L1:LINE\t5\nAluYb_001:AluYb:SINE\t7\n' > "$T/a.cntTable"
printf 'gene/TE\t/path/b.bam\nENSG0001\t0\nL1PA_0001:L1:LINE\t0\n' > "$T/b.cntTable"
mkdir -p "$T/custom/chimera/qc" "$T/custom/tecount/qc"
if ! python workflow/scripts/junction_qc_mqc.py \
      --tables "$T/jqc_a.tsv" "$T/jqc_b.tsv" --samples a b \
      --out "$T/custom/chimera/qc/junction_qc_mqc.json" \
      --out-te-gene-chimeras "$T/custom/chimera/qc/te_gene_chimeras_mqc.json" > "$T/jqc.log" 2>&1; then
  echo "ERROR: junction_qc_mqc.py failed"; cat "$T/jqc.log"; FAIL=1
fi
if ! python workflow/scripts/tecount_summary_mqc.py \
      --tables "$T/a.cntTable" "$T/b.cntTable" --samples a b \
      --out-assignment "$T/custom/tecount/qc/tecount_assignment_mqc.json" \
      --out-class "$T/custom/tecount/qc/tecount_te_class_mqc.json" > "$T/tsum.log" 2>&1; then
  echo "ERROR: tecount_summary_mqc.py failed"; cat "$T/tsum.log"; FAIL=1
fi
# A composition plot must offer exactly ONE counts/percentages control.
# cpswitch already provides it, and its percentage (share of the sample's
# total) is the right one -- also passing a hand-built percentage dataset in
# data_labels put a second, redundant switcher beside it.
python3 - "$T/custom/chimera/qc/junction_qc_mqc.json" <<'PY3' || FAIL=1
import json, sys
d = json.load(open(sys.argv[1]))
ok = True
if d.get("plot_type") == "bar":
    p = d.get("pconfig", {})
    if p.get("cpswitch") is not True:
        print("ERROR: composition plot should use cpswitch for its percentage view")
        ok = False
    if p.get("data_labels"):
        print(f"ERROR: cpswitch plus data_labels = two switchers; got "
              f"{[x.get('name') for x in p['data_labels']]}")
        ok = False
    if not isinstance(d.get("data"), dict):
        print("ERROR: a cpswitch plot takes ONE dataset, not a list")
        ok = False
sys.exit(0 if ok else 1)
PY3

if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/mqc" -n report "$T/custom" > "$T/mqc.log" 2>&1; then
  echo "ERROR: multiqc run failed"; tail -40 "$T/mqc.log"; FAIL=1
else
  for cid in chimera_junction_qc chimera_te_gene_chimeras \
             tecount_assignment tecount_te_class; do
    if ! grep -q "custom_content | $cid: Found" "$T/mqc.log" \
       || [ ! -f "$T/mqc/report_data/multiqc_${cid}_plot.txt" ]; then
      echo "ERROR: section '$cid' not rendered by MultiQC"
      grep "custom_content" "$T/mqc.log" | head -20
      ls "$T/mqc/report_data" | grep "${cid}" || true
      FAIL=1
    fi
  done
fi

exit $FAIL
