#!/usr/bin/env bash
# Guard 77: the SJ screen writes a per-pair count matrix keyed like the
# assembly screen's
#
# For differential analysis only the assembly screen had a matrix at the
# testable grain (gene_id:te_id:chimera_type); the SJ screen's matrix was
# per junction. aggregate_chimera_splice_junctions_counts.py sums junctions
# of one (gene, TE, type) per sample, keeps types apart, drops untyped or
# gene-less junctions, and writes an annotation table with the same key.
# The rule must be wired into the workflow and on by default.
#
# Run on its own:   .tests/guards/77_splice_junctions_gene_te_chimera_counts.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

S=workflow/scripts/aggregate_chimera_splice_junctions_counts.py
if [ ! -f "$S" ]; then
  echo "ERROR: $S is missing"; exit 1
fi

# te-gene-junctions.tsv.gz shape (only the columns the script reads)
{
  printf 'event_id\tgene_id\tte_id\tte_subfamily\tte_family\tte_class\tchimera_type\n'
  printf 'j1\tGENE1\tTE_A\tMTB\tERVL-MaLR\tLTR\tte_initiated\n'
  printf 'j2\tGENE1\tTE_A\tMTB\tERVL-MaLR\tLTR\tte_initiated\n'
  printf 'j3\tGENE1\tTE_A\tMTB\tERVL-MaLR\tLTR\tannotated_splice\n'
  printf 'j4\t.\tTE_B\tL1Md\tL1\tLINE\tte_initiated\n'
  printf 'j5\tGENE2\tTE_C\tB1\tAlu\tSINE\t.\n'
} | gzip -c > "$T/junctions.tsv.gz"
{
  printf 'event_id\tS1\tS2\n'
  printf 'j1\t3\t0\n'
  printf 'j2\t4\t6\n'
  printf 'j3\t100\t50\n'
  printf 'j4\t9\t9\n'
  printf 'j5\t7\t7\n'
} | gzip -c > "$T/counts.tsv.gz"
printf 'gene_id\tgene_name\nGENE1\tSpin1\n' | gzip -c > "$T/names.tsv.gz"
printf 'chr13\t100\t900\tGENE1\t.\t+\n' > "$T/genes.bed"
printf 'chr13\t200\t300\tTE_A\t.\t+\tERVL-MaLR\tLTR\tMTB\n' > "$T/te.bed"

if ! python3 "$S" --junctions "$T/junctions.tsv.gz" --counts "$T/counts.tsv.gz" \
      --gene-names "$T/names.tsv.gz" --genes-bed "$T/genes.bed" --te-bed "$T/te.bed" \
      --out-counts "$T/m.tsv.gz" --out-annotation "$T/a.tsv.gz" > "$T/log" 2>&1; then
  echo "ERROR: $S failed"; cat "$T/log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
def rows(p):
    with gzip.open(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

m = {r["gene_te_chimera_id"]: r for r in rows(f"{T}/m.tsv.gz")}
check(set(m) == {"GENE1:TE_A:te_initiated", "GENE1:TE_A:annotated_splice"},
      f"one row per typed (gene, TE, type), gene-less / untyped junctions "
      f"dropped; got {sorted(m)}")
ti = m.get("GENE1:TE_A:te_initiated", {})
check((ti.get("S1"), ti.get("S2")) == ("7", "6"),
      f"junctions of one (gene, TE, type) are summed per sample (want 7/6); got {ti}")
an = m.get("GENE1:TE_A:annotated_splice", {})
check((an.get("S1"), an.get("S2")) == ("100", "50"),
      f"a second type of the same pair stays its own row; got {an}")

a = {r["gene_te_chimera_id"]: r for r in rows(f"{T}/a.tsv.gz")}
check(set(a) == set(m), f"annotation keys must match counts keys; got {sorted(a)}")
r = a.get("GENE1:TE_A:te_initiated", {})
check(r.get("gene_symbol") == "Spin1" and r.get("gene_locus") == "chr13:101-900"
      and r.get("te_locus") == "chr13:201-300" and r.get("te_class") == "LTR"
      and r.get("n_junctions") == "2" and r.get("sj_event_ids") == "j1,j2",
      f"annotation row must carry symbol, 1-based loci, TE class and the "
      f"member junctions; got {r}")
sys.exit(0 if ok else 1)
PY

# wiring: the rule exists, declares the module it imports, and is on by default
R=workflow/rules/chimera_splice_junctions.smk
grep -q "rule chimera_splice_junctions_aggregate_counts" "$R" \
  || { echo "ERROR: no chimera_splice_junctions_aggregate_counts rule in $R"; FAIL=1; }
grep -q "aggregate_chimera_assembly_counts.py" "$R" \
  || { echo "ERROR: the rule must declare aggregate_chimera_assembly_counts.py (imported module)"; FAIL=1; }
python3 - <<'PY' || FAIL=1
import sys, yaml
d = yaml.safe_load(open("workflow/default-config/chimera.yaml"))
v = d["chimera"]["splice_junctions"]["outputs"].get("write_gene_te_chimera_counts")
if v is not True:
    print(f"ERROR: chimera.splice_junctions.outputs.write_gene_te_chimera_counts must default to true; got {v!r}")
    sys.exit(1)
PY

exit $FAIL
