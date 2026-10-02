#!/usr/bin/env bash
# Guard 78: the assembly screen names the TE the transcript actually uses
#
# Several TEs can hit one exon (hits are found with +/- breakpoint
# tolerance), and the classifier named the one with the lowest coordinate.
# On a real 84-sample run that named, e.g., an L1 ENDING at an MT2_Mm LTR
# promoter's TSS instead of the MT2_Mm the first exon runs through. Now:
#   - first exon: the TE containing the TSS, then most bases in the exon;
#   - internal exon: the TE with most bases in the exon;
#   - last exon: the TE containing the transcript's 3' end; with no TE
#     there, no single TE terminates it, so the old (coordinate) pick stays.
#
# Run on its own:   .tests/guards/78_assembly_names_the_te_at_the_transcript_end.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# GENE1 (+): exons [1000,1200) [2000,2200) [3000,3500)
printf 'chr1\t1000\t3500\tGENE1\t.\t+\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t3000\t3500\tGENE1\t.\t+\n' > "$T/exons.bed"
printf 'chr1\t3000\t3500\tTX1\t.\t+\n' > "$T/last_exons.bed"
: > "$T/first_exons.bed"
{
  # A_init: first exon [500,700), TSS 500. TE_UPX ends at the TSS (0 bp in
  # the exon), TE_PROM starts there and runs through it.
  printf 'chr1\t300\t500\tTE_UPX\t.\t+\tL1\tLINE\tL1Md\n'
  printf 'chr1\t500\t650\tTE_PROM\t.\t+\tERVL-MaLR\tLTR\tMT2_Mm\n'
  # A_exon: internal exon [1500,1600). TE_EDGE reaches 3 bp into it,
  # TE_EX lies inside it.
  printf 'chr1\t1400\t1503\tTE_EDGE\t.\t+\tAlu\tSINE\tB1\n'
  printf 'chr1\t1520\t1590\tTE_EX\t.\t+\tB2\tSINE\tB2_Mm1a\n'
  # A_term: last exon [4000,4800), 3' end 4800. TE_MID mid-UTR, TE_END
  # contains the 3' end.
  printf 'chr1\t4100\t4300\tTE_MID\t.\t+\tAlu\tSINE\tB1\n'
  printf 'chr1\t4700\t4900\tTE_END\t.\t+\tB2\tSINE\tB2_Mm2\n'
  # A_utr: last exon [5000,5800), no TE at the 3' end: TE_SMALL (lower
  # coordinate) keeps the name over the larger TE_BIG.
  printf 'chr1\t5100\t5200\tTE_SMALL\t.\t+\tAlu\tSINE\tB1\n'
  printf 'chr1\t5300\t5600\tTE_BIG\t.\t+\tL1\tLINE\tL1Md\n'
} > "$T/te.bed"

{
  printf 'chr1\tS\ttranscript\t501\t2200\t.\t+\t.\ttranscript_id "A_init"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t501\t700\t.\t+\t.\ttranscript_id "A_init"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t2001\t2200\t.\t+\t.\ttranscript_id "A_init"; gene_id "M.1";\n'
  printf 'chr1\tS\ttranscript\t1001\t2200\t.\t+\t.\ttranscript_id "A_exon"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t1001\t1200\t.\t+\t.\ttranscript_id "A_exon"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t1501\t1600\t.\t+\t.\ttranscript_id "A_exon"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t2001\t2200\t.\t+\t.\ttranscript_id "A_exon"; gene_id "M.2";\n'
  printf 'chr1\tS\ttranscript\t1001\t4800\t.\t+\t.\ttranscript_id "A_term"; gene_id "M.3";\n'
  printf 'chr1\tS\texon\t1001\t1200\t.\t+\t.\ttranscript_id "A_term"; gene_id "M.3";\n'
  printf 'chr1\tS\texon\t4001\t4800\t.\t+\t.\ttranscript_id "A_term"; gene_id "M.3";\n'
  printf 'chr1\tS\ttranscript\t1001\t5800\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t1001\t1200\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t5001\t5800\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.4";\n'
} > "$T/stringtie.gtf"

if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --last-exons "$T/last_exons.bed" \
      --te "$T/te.bed" --require-tss-in-te --breakpoint-tolerance 5 \
      --out "$T/asm.tsv" > "$T/asm.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/asm.log"; exit 1
fi

python3 - "$T" <<'PY' || FAIL=1
import csv, sys
T = sys.argv[1]
ok = True
def check(c, m):
    global ok
    if not c:
        print("ERROR:", m); ok = False
with open(f"{T}/asm.tsv") as fh:
    asm = {r["transcript_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
want = {
    "A_init": ("te_initiated", "TE_PROM",
               "first exon: the TE the exon runs through, not one ending at the TSS"),
    "A_exon": ("te_exonized", "TE_EX",
               "internal exon: the TE with most bases in the exon"),
    "A_term": ("te_terminated", "TE_END",
               "last exon: the TE containing the transcript's 3' end"),
    "A_utr": ("te_terminated", "TE_SMALL",
              "last exon with no TE at the 3' end keeps the coordinate pick"),
}
for tid, (ctype, te_id, why) in want.items():
    r = asm.get(tid, {})
    got = (r.get("chimera_type"), r.get("te_id"))
    check(got == (ctype, te_id), f"{tid}: {why}; want {(ctype, te_id)}, got {got}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
