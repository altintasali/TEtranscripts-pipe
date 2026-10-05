#!/usr/bin/env bash
# Guard 60: chimera_splice_junctions_summary_mqc.py's own output, and its
# rule appearing/disappearing with chimera.splice_junctions.enabled.
#
# Before this script existed, the splice_junctions screen had NO always-on
# report presence at all (its only two views, splice_junctions_pca_*/
# heatmap_*, are gated behind chimera.splice_junctions.qc.enabled, default
# false) -- a run with the screen enabled but that flag left at its default
# showed nothing for it.
#
# Run on its own:   .tests/guards/60_chimera_splice_junctions_summary_mqc_py_and_wiring.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- part 1: script-level output ---------------------------------------
mkdir -p "$T/sjs/per_sample"
cohort_header="event_id\tchrom\tintron_start\tintron_end\tstrand\tmotif\tcanonical\tannotated\tmulti_reads\toverhang\tdonor_hits\tacceptor_hits\tdirection\tdirection_ambiguous\tgene_id\tgene_strand\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\tte_initiated_detail\tantisense_flag\tlibrary_strand\ttranscript_strand\tgene_strand_match\tn_samples\ttotal_reads"
printf "$cohort_header\n" > "$T/sjs/cohort.tsv"
# E1: te_to_gene, te_initiated/upstream. E2: te_to_gene, te_initiated/internal
# (the exact classification-fix shape). E3: gene_to_te, te_exonized.
printf "E1\tchr1\t1\t2\t+\t1\tyes\t1\t0\t20\t.\t.\tte_to_gene\tno\tG1\t+\tT1\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\t2\t10\n" >> "$T/sjs/cohort.tsv"
printf "E2\tchr1\t3\t4\t+\t1\tyes\t1\t0\t20\t.\t.\tte_to_gene\tno\tG2\t+\tT2\t.\t.\t.\tte_initiated\tinternal\t.\tNA\tNA\tNA\t1\t5\n" >> "$T/sjs/cohort.tsv"
printf "E3\tchr1\t5\t6\t+\t1\tyes\t1\t0\t20\t.\t.\tgene_to_te\tno\tG3\t+\tT3\t.\t.\t.\tte_exonized\t.\t.\tNA\tNA\tNA\t1\t3\n" >> "$T/sjs/cohort.tsv"

per_header="event_id\tsample\tchrom\tintron_start\tintron_end\tstrand\tmotif\tcanonical\tannotated\tunique_reads\tmulti_reads\toverhang\tdonor_hits\tacceptor_hits\tdirection\tdirection_ambiguous\tgene_id\tgene_strand\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\tte_initiated_detail\tantisense_flag\tlibrary_strand\ttranscript_strand\tgene_strand_match"
printf "$per_header\n" > "$T/sjs/per_sample/S1_junctions_te-gene-junctions.tsv"
printf "E1\tS1\tchr1\t1\t2\t+\t1\tyes\t1\t7\t0\t20\t.\t.\tte_to_gene\tno\tG1\t+\tT1\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\n" >> "$T/sjs/per_sample/S1_junctions_te-gene-junctions.tsv"
printf "$per_header\n" > "$T/sjs/per_sample/S2_junctions_te-gene-junctions.tsv"
printf "E1\tS2\tchr1\t1\t2\t+\t1\tyes\t1\t3\t0\t20\t.\t.\tte_to_gene\tno\tG1\t+\tT1\t.\t.\t.\tte_initiated\tupstream\t.\tNA\tNA\tNA\n" >> "$T/sjs/per_sample/S2_junctions_te-gene-junctions.tsv"
printf "E2\tS2\tchr1\t3\t4\t+\t1\tyes\t1\t5\t0\t20\t.\t.\tte_to_gene\tno\tG2\t+\tT2\t.\t.\t.\tte_initiated\tinternal\t.\tNA\tNA\tNA\n" >> "$T/sjs/per_sample/S2_junctions_te-gene-junctions.tsv"
printf "E3\tS2\tchr1\t5\t6\t+\t1\tyes\t1\t3\t0\t20\t.\t.\tgene_to_te\tno\tG3\t+\tT3\t.\t.\t.\tte_exonized\t.\t.\tNA\tNA\tNA\n" >> "$T/sjs/per_sample/S2_junctions_te-gene-junctions.tsv"

if ! python3 workflow/scripts/chimera_splice_junctions_summary_mqc.py \
      --te-events "$T/sjs/cohort.tsv" \
      --per-sample "$T/sjs/per_sample/S1_junctions_te-gene-junctions.tsv" \
                   "$T/sjs/per_sample/S2_junctions_te-gene-junctions.tsv" \
      --sample-names S1 S2 \
      --out-highlights "$T/sjs/highlights.json" \
      --out-te-type "$T/sjs/te_type.json" > "$T/sjs/log" 2>&1; then
  echo "ERROR: chimera_splice_junctions_summary_mqc.py failed"; cat "$T/sjs/log"; FAIL=1
else
  python3 - "$T/sjs" <<'PY' || FAIL=1
import json, sys
d = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False

h = json.load(open(f"{d}/highlights.json"))
t = json.load(open(f"{d}/te_type.json"))
check(h.get("parent_id") == "chimera", f"highlights parent_id wrong: {h.get('parent_id')!r}")
check(t.get("parent_id") == "chimera", f"te_type parent_id wrong: {t.get('parent_id')!r}")
check(h["section_name"] == "Splice junctions - what this screen sees",
      f"highlights section misnamed: {h['section_name']!r}")
check(t["section_name"] == "Splice junctions - TE type",
      f"te_type section misnamed: {t['section_name']!r}")
check("100%" in h["description"] and "require_canonical" in h["description"],
      "highlights must say canonical rate is 100% by construction (require_canonical default)")
# per-sample, NOT deduplicated: E1 appears in both S1 and S2 -> counted in both
check(h["data"]["S1"]["te_to_gene"] == 1, f"S1 te_to_gene count wrong: {h['data']['S1']}")
check(h["data"]["S2"]["te_to_gene"] == 2, f"S2 te_to_gene count wrong: {h['data']['S2']}")
check(h["data"]["S2"]["gene_to_te"] == 1, f"S2 gene_to_te count wrong: {h['data']['S2']}")
# cohort-level, deduplicated: 3 unique events total
check("3" in h["description"], f"cohort total (3 unique events) not quoted: {h['description']}")
# THE POINT: te_initiated split into upstream/internal, not lumped together
te_data = t["data"]
check("TE-initiated (upstream)" in te_data and te_data["TE-initiated (upstream)"]["count"] == 1,
      f"te_initiated upstream count wrong: {te_data}")
check("TE-initiated (internal, skips an earlier exon)" in te_data
      and te_data["TE-initiated (internal, skips an earlier exon)"]["count"] == 1,
      f"te_initiated internal count wrong: {te_data}")
check("TE-exonized" in te_data and te_data["TE-exonized"]["count"] == 1,
      f"te_exonized count wrong: {te_data}")
sys.exit(0 if ok else 1)
PY
fi

# --- part 2: rule appears/disappears with chimera.splice_junctions.enabled
# config/test.yaml already has splice_junctions.enabled: true.
if ! snakemake --configfile config/test.yaml -n --cores 1 > "$T/on.log" 2>&1; then
  echo "ERROR: dry run with splice_junctions enabled failed"; tail -30 "$T/on.log"; FAIL=1
elif ! grep -qE "^chimera_splice_junctions_summary_mqc\b" "$T/on.log"; then
  echo "ERROR: chimera_splice_junctions_summary_mqc not scheduled with the screen enabled"
  FAIL=1
fi

sed '/^  splice_junctions:/,/^    enabled: true$/ s/^    enabled: true$/    enabled: false/' \
    config/test.yaml > "$T/off.yaml"
if ! grep -qE "^    enabled: false" "$T/off.yaml"; then
  echo "ERROR: guard 60 sed did not turn splice_junctions off"
  sed -n '/^  splice_junctions:/,+2p' "$T/off.yaml"; FAIL=1
elif ! snakemake --configfile "$T/off.yaml" -n --cores 1 > "$T/off.log" 2>&1; then
  echo "ERROR: dry run with splice_junctions disabled failed"; tail -30 "$T/off.log"; FAIL=1
elif grep -qE "^chimera_splice_junctions_summary_mqc\b" "$T/off.log"; then
  echo "ERROR: chimera_splice_junctions_summary_mqc still scheduled with the screen disabled"
  FAIL=1
fi

exit $FAIL
