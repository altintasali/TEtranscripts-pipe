#!/usr/bin/env bash
# Guard 78: the assembly screen names the TE the transcript actually uses
#
# Several TEs can hit one exon (hits are found with +/- breakpoint
# tolerance), and the classifier named the one with the lowest coordinate.
# On a real 84-sample run that named, e.g., an L1 ENDING at an MT2_Mm LTR
# promoter's TSS instead of the MT2_Mm the first exon runs through. Now:
#   - first exon: the TE containing the TSS, then most bases in the exon;
#   - internal exon: the TE with most bases in the exon;
#   - last exon: the TE nearest the last exon's splice acceptor -- the gene
#     splicing into a TE-derived last exon is what makes it TE-terminated,
#     and it is the TE the junction-based screens see. (Naming the TE at
#     the transcript's 3' end instead cut other-screen confirmation of the
#     renamed terminal calls to 7-15%, vs 41-46% for this pick.)
#
# Run on its own:   .tests/guards/78_assembly_names_the_te_the_transcript_uses.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# GENE1 (+): exons [1000,1200) [2000,2200) [3000,3500)
# GENE2 (-): exons [8000,8200) [9000,9200)
printf 'chr1\t1000\t3500\tGENE1\t.\t+\nchr1\t8000\t9200\tGENE2\t.\t-\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t3000\t3500\tGENE1\t.\t+\nchr1\t8000\t8200\tGENE2\t.\t-\nchr1\t9000\t9200\tGENE2\t.\t-\n' > "$T/exons.bed"
printf 'chr1\t3000\t3500\tTX1\t.\t+\nchr1\t8000\t8200\tTX2\t.\t-\n' > "$T/last_exons.bed"
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
  # A_term (+): last exon [4000,4800), acceptor 4000. TE_SPLIT ends 2 bp
  # into the exon, TE_ACC starts at the acceptor and runs into it, TE_END
  # holds the 3' end.
  printf 'chr1\t3800\t4002\tTE_SPLIT\t.\t+\tL2\tLINE\tL2a\n'
  printf 'chr1\t4000\t4200\tTE_ACC\t.\t+\tB2\tSINE\tB2_Mm2\n'
  printf 'chr1\t4700\t4900\tTE_END\t.\t+\tAlu\tSINE\tB1\n'
  # A_minus (-): last exon [7000,7800), acceptor 7800 (its 5' boundary on
  # "-"). TE_LOW sits at the 3' end (lowest coordinate), TE_HIGH 50 bp
  # from the acceptor.
  printf 'chr1\t7000\t7100\tTE_LOW\t.\t-\tAlu\tSINE\tB1\n'
  printf 'chr1\t7600\t7750\tTE_HIGH\t.\t-\tERVK\tLTR\tRLTR10\n'
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
  printf 'chr1\tS\ttranscript\t7001\t9200\t.\t-\t.\ttranscript_id "A_minus"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t7001\t7800\t.\t-\t.\ttranscript_id "A_minus"; gene_id "M.4";\n'
  printf 'chr1\tS\texon\t9001\t9200\t.\t-\t.\ttranscript_id "A_minus"; gene_id "M.4";\n'
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
    "A_term": ("te_terminated", "TE_ACC",
               "last exon: the TE at the splice acceptor running into the exon, "
               "not one ending there or the one at the 3' end"),
    "A_minus": ("te_terminated", "TE_HIGH",
                "last exon on '-': the TE nearest the acceptor, not the "
                "lowest-coordinate one at the 3' end"),
}
for tid, (ctype, te_id, why) in want.items():
    r = asm.get(tid, {})
    got = (r.get("chimera_type"), r.get("te_id"))
    check(got == (ctype, te_id), f"{tid}: {why}; want {(ctype, te_id)}, got {got}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
