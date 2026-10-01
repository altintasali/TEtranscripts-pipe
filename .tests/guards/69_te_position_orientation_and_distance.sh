#!/usr/bin/env bash
# Guard 69: candidates carry TE position, distance and orientation
#
# chimera_evidence.py adds, per gene-TE pair and from the annotation alone:
#   te_position          upstream / intronic / exonic / downstream, strand-aware
#   te_gene_distance_bp  gap to the gene span, 0 inside it
#   te_orientation       sense / antisense (TE strand vs gene strand)
# and the Candidates table shows them as a "TE vs gene" block (TE position,
# TE orientation, Distance) plus a Types column: the union of the chimera
# calls any screen made for the pair. Checked on both gene strands and all
# four positions, since "upstream" flips with the gene's strand.
#
# Run on its own:   .tests/guards/69_te_position_orientation_and_distance.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# GENE_P (+) [10000,20000), GENE_M (-) [50000,60000); three exons each.
printf 'chr1\t10000\t20000\tGENE_P\t.\t+\nchr1\t50000\t60000\tGENE_M\t.\t-\n' > "$T/genes.bed"
{
  printf 'chr1\t10000\t11000\tGENE_P\t.\t+\nchr1\t15000\t16000\tGENE_P\t.\t+\nchr1\t19000\t20000\tGENE_P\t.\t+\n'
  printf 'chr1\t50000\t51000\tGENE_M\t.\t-\nchr1\t55000\t56000\tGENE_M\t.\t-\nchr1\t59000\t60000\tGENE_M\t.\t-\n'
} > "$T/exons.bed"
# name -> expected (position, distance, orientation)
#   + gene: upstream = lower coordinates
#   TE_PU [8000,9000)   +  upstream    1000 sense
#   TE_PI [12000,12500) -  intronic    0    antisense
#   TE_PE [15500,15700) +  exonic      0    sense
#   TE_PD [23000,23300) -  downstream  3000 antisense
#   - gene: upstream = higher coordinates
#   TE_MU [62000,62500) -  upstream    2000 sense
#   TE_MI [52000,52300) -  intronic    0    sense
#   TE_ME [55100,55200) +  exonic      0    antisense
#   TE_MD [45000,46000) +  downstream  4000 antisense
{
  printf 'chr1\t8000\t9000\tTE_PU\t.\t+\tERVK\tLTR\tRLTR10\n'
  printf 'chr1\t12000\t12500\tTE_PI\t.\t-\tL1\tLINE\tL1Md\n'
  printf 'chr1\t15500\t15700\tTE_PE\t.\t+\tB2\tSINE\tB2_Mm1\n'
  printf 'chr1\t23000\t23300\tTE_PD\t.\t-\tERVL\tLTR\tMT2A\n'
  printf 'chr1\t62000\t62500\tTE_MU\t.\t-\tERVL-MaLR\tLTR\tMTA\n'
  printf 'chr1\t52000\t52300\tTE_MI\t.\t-\tL1\tLINE\tL1Md\n'
  printf 'chr1\t55100\t55200\tTE_ME\t.\t+\tAlu\tSINE\tB1\n'
  printf 'chr1\t45000\t46000\tTE_MD\t.\t+\tERVK\tLTR\tIAP\n'
} > "$T/te.bed"

# Merged SJ table: one te_initiated call per pair; GENE_P/TE_PU also has an
# annotated splice, which Types must leave out (not a chimera call).
python3 - "$T" <<'PY'
import gzip, sys
T = sys.argv[1]
cols = ["event_id", "gene_id", "te_id", "te_subfamily", "te_family", "te_class",
        "canonical", "chimera_type", "n_samples", "total_reads"]
pairs = [("GENE_P", t) for t in ("TE_PU", "TE_PI", "TE_PE", "TE_PD")] + \
        [("GENE_M", t) for t in ("TE_MU", "TE_MI", "TE_ME", "TE_MD")]
with gzip.open(f"{T}/sj.tsv.gz", "wt") as fh:
    fh.write("\t".join(cols) + "\n")
    for i, (g, t) in enumerate(pairs):
        fh.write("\t".join([f"s{i}", g, t, "sf", "fam", "LTR", "yes",
                            "te_initiated", "2", "5"]) + "\n")
    fh.write("\t".join(["s_annot", "GENE_P", "TE_PU", "sf", "fam", "LTR", "yes",
                        "annotated_splice", "4", "500"]) + "\n")
with gzip.open(f"{T}/cr.tsv.gz", "wt") as fh:
    fh.write("event_id\tgene_id\tte_id\tchimera_type\n")
PY

if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --sj "$T/sj.tsv.gz" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --out "$T/cand.tsv.gz" > "$T/ev.log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev.log"; exit 1
fi
# without the annotation the three columns degrade to "." instead of failing
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --sj "$T/sj.tsv.gz" --out "$T/cand_noann.tsv.gz" > "$T/ev2.log" 2>&1; then
  echo "ERROR: chimera_evidence.py without --genes/--exons/--te failed"; cat "$T/ev2.log"; exit 1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/cand.tsv.gz" --out "$T/table_mqc.json" > "$T/tbl.log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed"; cat "$T/tbl.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, json, sys
T = sys.argv[1]
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
def rows(path):
    with gzip.open(path, "rt") as fh:
        return {(r["gene_id"], r["te_id"]): r for r in csv.DictReader(fh, delimiter="\t")}

want = {
    ("GENE_P", "TE_PU"): ("upstream", "1000", "sense"),
    ("GENE_P", "TE_PI"): ("intronic", "0", "antisense"),
    ("GENE_P", "TE_PE"): ("exonic", "0", "sense"),
    ("GENE_P", "TE_PD"): ("downstream", "3000", "antisense"),
    ("GENE_M", "TE_MU"): ("upstream", "2000", "sense"),
    ("GENE_M", "TE_MI"): ("intronic", "0", "sense"),
    ("GENE_M", "TE_ME"): ("exonic", "0", "antisense"),
    ("GENE_M", "TE_MD"): ("downstream", "4000", "antisense"),
}
cand = rows(f"{T}/cand.tsv.gz")
for pair, (pos, dist, ori) in want.items():
    r = cand.get(pair, {})
    got = (r.get("te_position"), r.get("te_gene_distance_bp"), r.get("te_orientation"))
    check(got == (pos, dist, ori),
          f"{pair}: want position/distance/orientation {(pos, dist, ori)}, got {got}")
some = next(iter(cand.values()))
check("cr_gene_te_distance" not in some,
      "cr_gene_te_distance was folded into te_gene_distance_bp and must be gone")

for pair, r in rows(f"{T}/cand_noann.tsv.gz").items():
    got = (r.get("te_position"), r.get("te_gene_distance_bp"), r.get("te_orientation"))
    check(got == (".", ".", "."), f"{pair} without annotation: want '.' x3, got {got}")

d = json.load(open(f"{T}/table_mqc.json"))
headers = d.get("headers", {})
for h in ("TE position", "TE orientation", "Distance", "Types"):
    check(h in headers, f"Candidates table is missing the {h!r} column")
# the TE-vs-gene block sits right after TE class, before Screens
place = {h: v.get("placement") for h, v in headers.items()}
check(place["TE class"] < place["TE position"] < place["TE orientation"]
      < place["Distance"] < place["Screens"],
      f"TE vs gene block must sit between TE class and Screens; placements {place}")
e = d["data"].get("GENE_P | TE_PU", {})
check(e.get("TE position") == "upstream" and e.get("TE orientation") == "sense"
      and e.get("Distance") == 1000,
      f"table must carry the pair's position/orientation/distance; got {e}")
check(e.get("Types") == "te_initiated",
      f"Types is the union of chimera CALLS only (annotated_splice left out); "
      f"got {e.get('Types')!r}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
