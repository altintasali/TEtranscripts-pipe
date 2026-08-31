#!/usr/bin/env bash
# Guard 46: each screen gets its own TE-type section, and both say the labels differ
#
# Run on its own:   .tests/guards/46_te_type_shows_both_screens_and_says_they_differ.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- Both screens emit te_initiated / te_exonized / te_terminated, and they
# are NOT the same measurement: the read screen classifies by where the TE
# sits relative to the gene body, the assembly screen by which exon of the
# assembled transcript it hits. They are therefore kept in SEPARATE sections
# -- one shared plot would invite reading agreement as corroboration -- and
# each must spell the difference out.
mkdir -p "$T/tt/qc"
for s in s1 s2; do
  printf 'metric\tvalue\nchimera_type_te_initiated\t10\nchimera_type_te_exonized\t4\nchimera_type_te_terminated\t2\n' \
    | gzip -c > "$T/tt/${s}_qc.tsv.gz"
done
{ printf 'transcript_id\tte_id\tmatched_gene_id\tstrand_match\tchimera_type\n'
  printf 'T1\tTE_A\tG1\tyes\tte_initiated\n'
  printf 'T2\tTE_B\tG2\tno\tte_exonized\n'
  printf 'T3\tTE_C\t.\tyes\tte_initiated_intergenic\n'
} | gzip -c > "$T/tt/cands.tsv.gz"

if ! python3 workflow/scripts/chimera_te_type_mqc.py \
      --qc-tables "$T/tt/s1_qc.tsv.gz" "$T/tt/s2_qc.tsv.gz" --samples s1 s2 \
      --out "$T/tt/qc/chimera_reads_te_type_mqc.json" > "$T/tt/log" 2>&1; then
  echo "ERROR: chimera_te_type_mqc.py failed"; cat "$T/tt/log"; FAIL=1
elif ! python3 workflow/scripts/chimera_assembly_summary_mqc.py \
      --candidates "$T/tt/cands.tsv.gz" \
      --out-classes "$T/tt/qc/chimera_assembly_classes_mqc.json" \
      --out-highlights "$T/tt/qc/chimera_assembly_highlights_mqc.json" \
      --out-strand-rate "$T/tt/qc/chimera_assembly_strand_rate_mqc.json" \
      >> "$T/tt/log" 2>&1; then
  echo "ERROR: chimera_assembly_summary_mqc.py failed"; cat "$T/tt/log"; FAIL=1
else
  python3 - "$T/tt/qc" <<'ASSERTEOF' || FAIL=1
import json, sys
d = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False

reads = json.load(open(f"{d}/chimera_reads_te_type_mqc.json"))
asm = json.load(open(f"{d}/chimera_assembly_classes_mqc.json"))

# one section per screen, both inside the single Chimera group
check(reads["section_name"] == "Reads - TE type",
      f"read section misnamed: {reads['section_name']!r}")
check(asm["section_name"] == "Assembly - TE type",
      f"assembly section misnamed: {asm['section_name']!r}")
for doc, name in ((reads, "reads"), (asm, "assembly")):
    check(doc.get("parent_id") == "chimera", f"{name} parent_id wrong")

# the read view is PER SAMPLE -- that view existed nowhere else once the
# per-sample status grid was removed
check(set(reads["data"]) == {"s1", "s2"},
      f"read TE type must be per sample; got {list(reads['data'])}")
check(reads["data"]["s1"]["te_initiated"] == 10,
      f"per-sample counts wrong: {reads['data']['s1']}")

# THE POINT: each section must say the other screen means something else by
# the same words, or a reader reads agreement as corroboration
for doc, name, other in ((reads, "Reads", "transcript structure"),
                         (asm, "Assembly", "genomic position")):
    desc = doc["description"]
    check("same words for a different measurement" in desc,
          f"{name} section must warn the labels are not comparable")
    check(other in desc,
          f"{name} section must describe the other screen's basis ({other})")
    check("neither is wrong" in desc.lower(),
          f"{name} section must say neither classification is wrong")
sys.exit(0 if ok else 1)
ASSERTEOF
fi

if ! multiqc --force --no-ansi -c workflow/default-config/multiqc_config.yaml \
      -o "$T/tt/out" -n r "$T/tt/qc" > "$T/tt/render.log" 2>&1; then
  echo "ERROR: multiqc failed on the TE type sections"
  tail -30 "$T/tt/render.log"; FAIL=1
else
  for want in "Reads - TE type" "Assembly - TE type"; do
    grep -q "$want" "$T/tt/out/r.html" || { echo "ERROR: '$want' not rendered"; FAIL=1; }
  done
fi

# an empty run must not take the whole report down
printf 'metric\tvalue\nchimera_type_te_initiated\t0\n' | gzip -c > "$T/tt/empty_qc.tsv.gz"
if ! python3 workflow/scripts/chimera_te_type_mqc.py \
      --qc-tables "$T/tt/empty_qc.tsv.gz" --samples s1 --out "$T/tt/empty_mqc.json" \
      > "$T/tt/empty.log" 2>&1; then
  echo "ERROR: chimera_te_type_mqc.py died on an empty run"; cat "$T/tt/empty.log"; FAIL=1
elif ! python3 -c "
import json,sys
sys.exit(0 if json.load(open('$T/tt/empty_mqc.json'))['plot_type']=='html' else 1)"; then
  echo "ERROR: an all-zero TE type plot must fall back to html"; FAIL=1
fi

exit $FAIL
