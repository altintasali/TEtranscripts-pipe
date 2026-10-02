#!/usr/bin/env bash
# Guard 67: all three chimera classifiers apply the strand rule
#
# A gene<->TE call whose transcript runs on the strand OPPOSITE its gene is
# antisense transcription through that gene's exon, not initiation /
# termination / exonization of the gene. On a real 4-sample run these were
# the te_initiated calls whose TE sat DOWNSTREAM of the whole gene (and the
# te_terminated ones with the TE upstream) -- geometrically impossible for a
# sense chimera. All three classifiers now:
#   - prefer a gene on the transcript's own strand when a breakpoint/exon
#     overlaps genes on both strands (prefer_gene_on_strand), and
#   - type a call still on the opposite strand antisense_to_gene.
# Also pinned here:
#   - SJ.out.tab's strand is the TRANSCRIPT strand (from the intron motif)
#     for any library -- it used to be flipped for reverse (dUTP) libraries,
#     inverting gene_strand_match there;
#   - the chimeric-reads table's gene_te_distance column ("trans" / 0 / bp).
#
# Run on its own:   .tests/guards/67_classifiers_apply_the_strand_rule.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# ---------------------------------------------------------------- fixtures
# block 1 (chr1 0-5k):    GENE1 (+), exons [1000,1200) [2000,2200); TE_UP [500,700); TE_DN [2500,2700)
# block 2 (chr1 10-12k):  GENE_M (-) and GENE_P (+) share exon [11000,11200); TE_OV [10500,10700)
#                         -> sorted order puts GENE_M first; the + transcript must still get GENE_P
# chr2:                   TE_FAR [500,700) for the trans case
printf 'chr1\t1000\t2200\tGENE1\t.\t+\nchr1\t11000\t11200\tGENE_M\t.\t-\nchr1\t11000\t11200\tGENE_P\t.\t+\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t11000\t11200\tGENE_M\t.\t-\nchr1\t11000\t11200\tGENE_P\t.\t+\n' > "$T/exons.bed"
printf 'chr1\t500\t700\tTE_UP\t.\t+\tERVL-MaLR\tLTR\tMTA\nchr1\t2500\t2700\tTE_DN\t.\t-\tERVL-MaLR\tLTR\tMTB\nchr1\t10500\t10700\tTE_OV\t.\t+\tERVL\tLTR\tMT2A\nchr2\t500\t700\tTE_FAR\t.\t+\tL1\tLINE\tL1Md\n' > "$T/te.bed"
: > "$T/first_exons.bed"

# ---------------------------------------------------------------- SJ
# SJ.out.tab: chrom, intron_start, intron_end, strand(1=+,2=-), motif,
# annotated, unique, multi, overhang. Donor/acceptor follow the strand.
#   S_plus : TE_UP -> GENE1 exon1 on + (gene strand)        -> te_initiated
#   S_minus: same intron, strand -: donor in GENE1 exon1,
#            acceptor in TE_UP, transcript opposite GENE1   -> antisense_to_gene
#   S_ovl  : TE_OV -> shared exon on +; GENE_M sorts first   -> GENE_P, te_initiated
printf 'chr1\t701\t1000\t1\t1\t0\t5\t0\t20\n'   >  "$T/sj.tab"
printf 'chr1\t701\t1000\t2\t2\t0\t5\t0\t20\n'   >> "$T/sj.tab"
printf 'chr1\t10701\t11000\t1\t1\t0\t5\t0\t20\n' >> "$T/sj.tab"
for lib in no reverse; do
  if ! python3 workflow/scripts/classify_chimera_splice_junctions.py \
        --sj "$T/sj.tab" --genes "$T/genes.bed" --exons "$T/exons.bed" \
        --te "$T/te.bed" --sample S1 --library-strandedness "$lib" \
        --out "$T/sj_$lib.tsv" > "$T/sj_$lib.log" 2>&1; then
    echo "ERROR: classify_chimera_splice_junctions.py ($lib) failed"; cat "$T/sj_$lib.log"; FAIL=1
  fi
done

# ---------------------------------------------------------------- CR
# Chimeric.out.junction: donor_chrom, donor_bp (first base of the donor's
# intron), donor_strand, acceptor_chrom, acceptor_bp (last base of the
# acceptor's intron), acceptor_strand, junction_type, repeat_left, ...
#   C_local: TE_UP -> GENE1 exon1, both segments +          -> te_initiated (forward lib)
#                                                             antisense_to_gene (reverse lib)
#   C_trans: TE_FAR (chr2) -> GENE1 exon1                   -> "." , distance "trans"
printf 'chr1\t700\t+\tchr1\t1000\t+\t1\t0\t0\tr1\n'  >  "$T/chim.junction"
printf 'chr2\t700\t+\tchr1\t1000\t+\t1\t0\t0\tr2\n'  >> "$T/chim.junction"
for lib in forward reverse no; do
  if ! python3 workflow/scripts/classify_chimera_chimeric_reads.py \
        --junctions "$T/chim.junction" --genes "$T/genes.bed" --exons "$T/exons.bed" \
        --te "$T/te.bed" --sample S1 --library-strandedness "$lib" \
        --out "$T/cr_$lib.tsv" > "$T/cr_$lib.log" 2>&1; then
    echo "ERROR: classify_chimera_chimeric_reads.py ($lib) failed"; cat "$T/cr_$lib.log"; FAIL=1
  fi
done

# ---------------------------------------------------------------- assembly
#   A_anti: '-' transcript, TSS inside TE_DN (first exon 2501-2700 on '-'),
#           downstream exon on GENE1 exon1 (+)              -> antisense_to_gene
#   A_ovl : '+' transcript from TE_OV into the shared exon;
#           GENE_M (-) would be hit first                    -> GENE_P, te_initiated
{
  printf 'chr1\tStringTie\ttranscript\t1001\t2700\t.\t-\t.\ttranscript_id "A_anti"; gene_id "M.1";\n'
  printf 'chr1\tStringTie\texon\t1001\t1200\t.\t-\t.\ttranscript_id "A_anti"; gene_id "M.1";\n'
  printf 'chr1\tStringTie\texon\t2501\t2700\t.\t-\t.\ttranscript_id "A_anti"; gene_id "M.1";\n'
  printf 'chr1\tStringTie\ttranscript\t10501\t11200\t.\t+\t.\ttranscript_id "A_ovl"; gene_id "M.2";\n'
  printf 'chr1\tStringTie\texon\t10501\t10700\t.\t+\t.\ttranscript_id "A_ovl"; gene_id "M.2";\n'
  printf 'chr1\tStringTie\texon\t11001\t11200\t.\t+\t.\ttranscript_id "A_ovl"; gene_id "M.2";\n'
} > "$T/stringtie.gtf"
if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --te "$T/te.bed" \
      --out "$T/asm.tsv" > "$T/asm.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/asm.log"; FAIL=1
fi

[ "$FAIL" = 0 ] || exit $FAIL

python3 - "$T" <<'PY' || FAIL=1
import csv, sys
T = sys.argv[1]
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

# ---- SJ
for lib in ("no", "reverse"):
    r = rows(f"{T}/sj_{lib}.tsv")
    by = {(x["intron_start"], x["strand"]): x for x in r}
    plus, minus, ovl = by[("701", "+")], by[("701", "-")], by[("10701", "+")]
    check(plus["chimera_type"] == "te_initiated",
          f"SJ[{lib}] sense junction TE->gene must stay te_initiated, got {plus['chimera_type']!r}")
    check(plus["transcript_strand"] == "+" and plus["gene_strand_match"] == "yes",
          f"SJ[{lib}] SJ.out.tab strand is the transcript strand for ANY library "
          f"(it used to be flipped for reverse); got transcript_strand="
          f"{plus['transcript_strand']!r}, gene_strand_match={plus['gene_strand_match']!r}")
    check(minus["chimera_type"] == "antisense_to_gene",
          f"SJ[{lib}] junction on the gene's opposite strand must be antisense_to_gene, "
          f"got {minus['chimera_type']!r}")
    check(minus["te_initiated_detail"] == ".",
          f"SJ[{lib}] antisense_to_gene carries no initiated detail, got {minus['te_initiated_detail']!r}")
    check(ovl["gene_id"] == "GENE_P" and ovl["chimera_type"] == "te_initiated",
          f"SJ[{lib}] overlapping genes: the + junction must pick the + gene GENE_P "
          f"(not first-sorted GENE_M) and stay te_initiated; got "
          f"{ovl['gene_id']!r}/{ovl['chimera_type']!r}")

# ---- CR
def cr(lib):
    r = rows(f"{T}/cr_{lib}.tsv")
    local = [x for x in r if x["donor_chrom"] == "chr1"][0]
    trans = [x for x in r if x["donor_chrom"] == "chr2"][0]
    return local, trans
local, trans = cr("forward")
check(local["chimera_type"] == "te_initiated" and local["gene_strand_match"] == "yes",
      f"CR[forward] sense event must stay te_initiated/match=yes, got "
      f"{local['chimera_type']!r}/{local['gene_strand_match']!r}")
check(local["gene_te_distance"] == "300",
      f"CR gene_te_distance for TE_UP [500,700) vs GENE1 [1000,2200) must be 300, "
      f"got {local['gene_te_distance']!r}")
check(trans["gene_te_distance"] == "trans" and trans["chimera_type"] == ".",
      f"CR trans event: distance 'trans', type '.'; got "
      f"{trans['gene_te_distance']!r}/{trans['chimera_type']!r}")
local, _ = cr("reverse")
check(local["chimera_type"] == "antisense_to_gene",
      f"CR[reverse] same reads on a reverse library are the gene's opposite strand "
      f"-> antisense_to_gene, got {local['chimera_type']!r}")
local, _ = cr("no")
check(local["chimera_type"] == "te_initiated" and local["gene_strand_match"] == "NA",
      f"CR[unstranded] strand unknown -> type unchanged, match NA; got "
      f"{local['chimera_type']!r}/{local['gene_strand_match']!r}")

# ---- assembly
r = {x["transcript_id"]: x for x in rows(f"{T}/asm.tsv")}
check("A_anti" in r and r["A_anti"]["chimera_type"] == "antisense_to_gene",
      f"assembly '-' transcript onto a '+' gene must be antisense_to_gene, got "
      f"{r.get('A_anti', {}).get('chimera_type')!r}")
check("A_ovl" in r and r["A_ovl"]["matched_gene_id"] == "GENE_P"
      and r["A_ovl"]["chimera_type"] == "te_initiated"
      and r["A_ovl"]["strand_match"] == "yes",
      f"assembly overlapping genes: '+' transcript must match GENE_P and stay "
      f"te_initiated; got {r.get('A_ovl')!r}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
