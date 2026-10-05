#!/usr/bin/env bash
# Guard 79: the assembly screen reports how far a last-exon TE sits from
# the exon's splice acceptor
#
# A TE-terminated call only needs some TE in a new last exon, and last exons
# can be long 3' UTRs. On a real 84-sample run, TEs 6-200 bp past the
# acceptor were 63-68% sense to the transcript (a TE that can end it), while
# those beyond ~500 bp sat at background. The distance is now a signal, not
# a filter:
#   - classify_chimera_assembly.py: te_acceptor_distance_bp (0 = the TE takes
#     the splice; "." for a TE in any other exon), on either strand;
#   - chimera_evidence.py: assembly_te_acceptor_distance_bp, the smallest
#     over the pair's calls;
#   - the Candidates table: "Assembly last-exon TE distance".
#
# Run on its own:   .tests/guards/79_assembly_last_exon_te_distance.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# GENE1 (+): exons [1000,1200) [2000,2200); GENE2 (-): [8000,8200) [9000,9200)
printf 'chr1\t1000\t2200\tGENE1\t.\t+\nchr1\t8000\t9200\tGENE2\t.\t-\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t8000\t8200\tGENE2\t.\t-\nchr1\t9000\t9200\tGENE2\t.\t-\n' > "$T/exons.bed"
printf 'chr1\t2000\t2200\tTX1\t.\t+\nchr1\t8000\t8200\tTX2\t.\t-\n' > "$T/last_exons.bed"
: > "$T/first_exons.bed"
{
  # T_at (+): last exon [4000,4800), acceptor 4000; TE_AT takes the splice
  printf 'chr1\t3990\t4100\tTE_AT\t.\t+\tB2\tSINE\tB2_Mm2\n'
  # T_near (+): last exon [5000,5800), acceptor 5000; TE_NEAR 120 bp past it
  printf 'chr1\t5120\t5300\tTE_NEAR\t.\t+\tB2\tSINE\tB2_Mm1a\n'
  # T_minus (-): last exon [7000,7800), acceptor 7800; TE_FAR ends 700 bp before it
  printf 'chr1\t7000\t7100\tTE_FAR\t.\t-\tAlu\tSINE\tB1\n'
  # T_init (+): first-exon TE, distance must be "."
  printf 'chr1\t500\t700\tTE_PROM\t.\t+\tERVL-MaLR\tLTR\tMT2_Mm\n'
} > "$T/te.bed"
{
  printf 'chr1\tS\ttranscript\t1001\t4800\t.\t+\t.\ttranscript_id "T_at"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T_at"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t4001\t4800\t.\t+\t.\ttranscript_id "T_at"; gene_id "M.1";\n'
  printf 'chr1\tS\ttranscript\t1001\t5800\t.\t+\t.\ttranscript_id "T_near"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T_near"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t5001\t5800\t.\t+\t.\ttranscript_id "T_near"; gene_id "M.2";\n'
  printf 'chr1\tS\ttranscript\t7001\t9200\t.\t-\t.\ttranscript_id "T_minus"; gene_id "M.3";\n'
  printf 'chr1\tS\texon\t7001\t7800\t.\t-\t.\ttranscript_id "T_minus"; gene_id "M.3";\n'
  printf 'chr1\tS\texon\t9001\t9200\t.\t-\t.\ttranscript_id "T_minus"; gene_id "M.3";\n'
  printf 'chr1\tS\ttranscript\t501\t2200\t.\t+\t.\ttranscript_id "T_init"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t501\t700\t.\t+\t.\ttranscript_id "T_init"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t2001\t2200\t.\t+\t.\ttranscript_id "T_init"; gene_id "M.4";\n'
} > "$T/stringtie.gtf"

if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --last-exons "$T/last_exons.bed" \
      --te "$T/te.bed" --require-tss-in-te --breakpoint-tolerance 5 \
      --out "$T/asm.tsv.gz" > "$T/asm.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/asm.log"; exit 1
fi
# a second TE_NEAR call farther out, to check the pair keeps the smallest
python3 - "$T" <<'PY'
import csv, gzip, sys
T = sys.argv[1]
with gzip.open(f"{T}/asm.tsv.gz", "rt") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = list(rows[0])
near = next(r for r in rows if r["transcript_id"] == "T_near")
rows.append(dict(near, transcript_id="T_near2", te_acceptor_distance_bp="900"))
with gzip.open(f"{T}/asm2.tsv.gz", "wt") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(r[c] for c in cols) + "\n")
PY
printf 'event_id\tgene_id\tte_id\tchimera_type\tn_samples\ttotal_reads\n' | gzip -c > "$T/cr.tsv.gz"
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --assembly "$T/asm2.tsv.gz" \
      --out "$T/cand.tsv.gz" > "$T/ev.log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev.log"; exit 1
fi
if ! python3 workflow/scripts/chimera_candidates_table_mqc.py \
      --evidence "$T/cand.tsv.gz" --out "$T/table.json" > "$T/tbl.log" 2>&1; then
  echo "ERROR: chimera_candidates_table_mqc.py failed"; cat "$T/tbl.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, json, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
with gzip.open(f"{T}/asm.tsv.gz", "rt") as fh:
    asm = {r["transcript_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
want = {"T_at": ("te_terminated", "0"), "T_near": ("te_terminated", "120"),
        "T_minus": ("te_terminated", "700"), "T_init": ("te_initiated", ".")}
for tid, w in want.items():
    r = asm.get(tid, {})
    got = (r.get("chimera_type"), r.get("te_acceptor_distance_bp"))
    check(got == w, f"{tid}: want (type, te_acceptor_distance_bp) {w}, got {got}")

with gzip.open(f"{T}/cand.tsv.gz", "rt") as fh:
    cand = {(r["gene_id"], r["te_id"]): r for r in csv.DictReader(fh, delimiter="\t")}
check(cand.get(("GENE1", "TE_NEAR"), {}).get("assembly_te_acceptor_distance_bp") == "120",
      f"pair column is the smallest over its calls (120, not 900); got "
      f"{cand.get(('GENE1', 'TE_NEAR'), {}).get('assembly_te_acceptor_distance_bp')!r}")
check(cand.get(("GENE1", "TE_PROM"), {}).get("assembly_te_acceptor_distance_bp") == ".",
      "a pair with no last-exon call must be '.'")

d = json.load(open(f"{T}/table.json"))
check("Assembly last-exon TE distance" in d.get("headers", {}),
      "Candidates table is missing 'Assembly last-exon TE distance'")
rows = d.get("data", {})
e = next((v for k, v in rows.items() if "TE_NEAR" in k), None)
if e is not None:
    check(e.get("Assembly last-exon TE distance") == 120,
          f"table must carry the distance; got {e}")
sys.exit(0 if ok else 1)
PY

grep -q '"Assembly last-exon TE distance"' workflow/scripts/chimera_candidates_explorer.R \
  || { echo "ERROR: explorer has no 'Assembly last-exon TE distance' column"; FAIL=1; }

exit $FAIL
