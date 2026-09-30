#!/usr/bin/env bash
# Guard 68: known gene structure is its own class, and far/trans
# chimeric-read events do not make candidates
#
# Measured on a real 4-sample run after the strand rule:
#   - 25% of SJ gene-TE junctions were annotated GTF introns -- ordinary
#     splicing whose splice site falls inside a TE, not a new chimera;
#   - 43% of assembly te_terminated calls had their TE in an annotated last
#     exon -- an ordinary 3' UTR TE, not a new TE-terminated transcript;
#   - 97.5% of chimeric-read gene-TE events were trans or >200 kb from the
#     gene -- the random-partner pattern of template switching / chimeric
#     ligation.
# The first two are now typed annotated_splice /
# annotated_terminal_exon_embedded_te (rows kept, but no longer
# te_initiated/te_terminated/te_exonized); the third no longer counts toward
# candidates.tsv.gz (chimera.chimeric_reads.max_gene_te_distance).
#
# Run on its own:   .tests/guards/68_known_structure_and_far_chimeric_reads_are_not_candidates.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# ------------------------------------------------ reference: GTF -> features
# GENE1 (+): TX1 exons 1001-1200, 2001-2200, 3001-3500 (last exon has TE_UTR)
{
  printf 'chr1\tx\texon\t1001\t1200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n'
  printf 'chr1\tx\texon\t2001\t2200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n'
  printf 'chr1\tx\texon\t3001\t3500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n'
} > "$T/genes.gtf"
if ! python3 workflow/scripts/annotation_splice_features.py --gtf "$T/genes.gtf" \
      --out-introns "$T/introns.tsv.gz" --out-last-exons "$T/last_exons.bed" \
      > "$T/feat.log" 2>&1; then
  echo "ERROR: annotation_splice_features.py failed"; cat "$T/feat.log"; exit 1
fi

printf 'chr1\t1000\t3500\tGENE1\t.\t+\n' > "$T/genes.bed"
printf 'chr1\t1000\t1200\tGENE1\t.\t+\nchr1\t2000\t2200\tGENE1\t.\t+\nchr1\t3000\t3500\tGENE1\t.\t+\n' > "$T/exons.bed"
# TE_IN1 straddles exon1's 3' end (so annotated intron 1201-2000 has its
# donor inside a TE); TE_INTRON sits in intron 2; TE_UTR in the last exon.
printf 'chr1\t1150\t1250\tTE_IN1\t.\t+\tERVL\tLTR\tMT2A\nchr1\t2500\t2700\tTE_INTRON\t.\t+\tERVL-MaLR\tLTR\tMTA\nchr1\t3200\t3400\tTE_UTR\t.\t+\tAlu\tSINE\tB1\n' > "$T/te.bed"
: > "$T/first_exons.bed"

# ------------------------------------------------ SJ
#   J_annot: 1201-2000 (+) -- the annotated intron exon1->exon2; donor window
#            overlaps TE_IN1 and exon1 -> must be annotated_splice
#   J_novel: 2701-3000 (+) -- TE_INTRON -> exon3, not annotated -> te_initiated
printf 'chr1\t1201\t2000\t1\t1\t1\t9\t0\t30\nchr1\t2701\t3000\t1\t1\t0\t5\t0\t20\n' > "$T/sj.tab"
if ! python3 workflow/scripts/classify_chimera_splice_junctions.py \
      --sj "$T/sj.tab" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --te "$T/te.bed" --annotated-introns "$T/introns.tsv.gz" --sample S1 \
      --out "$T/sj.tsv" > "$T/sj.log" 2>&1; then
  echo "ERROR: classify_chimera_splice_junctions.py failed"; cat "$T/sj.log"; FAIL=1
fi

# ------------------------------------------------ assembly
#   A_utr  : exon2 -> last exon 3001-3500 containing TE_UTR (= annotated
#            last exon)                         -> annotated_terminal_exon_embedded_te
#   A_novel: exon2 -> novel last exon 2501-2700 inside TE_INTRON -> te_terminated
{
  printf 'chr1\tS\ttranscript\t2001\t3500\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t2001\t2200\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.1";\n'
  printf 'chr1\tS\texon\t3001\t3500\t.\t+\t.\ttranscript_id "A_utr"; gene_id "M.1";\n'
  printf 'chr1\tS\ttranscript\t2001\t2700\t.\t+\t.\ttranscript_id "A_novel"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t2001\t2200\t.\t+\t.\ttranscript_id "A_novel"; gene_id "M.2";\n'
  printf 'chr1\tS\texon\t2501\t2700\t.\t+\t.\ttranscript_id "A_novel"; gene_id "M.2";\n'
} > "$T/stringtie.gtf"
if ! python3 workflow/scripts/classify_chimera_assembly.py \
      --gtf "$T/stringtie.gtf" --genes "$T/genes.bed" --exons "$T/exons.bed" \
      --first-exons "$T/first_exons.bed" --last-exons "$T/last_exons.bed" \
      --te "$T/te.bed" --out "$T/asm.tsv" > "$T/asm.log" 2>&1; then
  echo "ERROR: classify_chimera_assembly.py failed"; cat "$T/asm.log"; FAIL=1
fi

# ------------------------------------------------ candidates: far CR events
# merged CR table (chimera_chimeric_reads_counts.py shape, relevant columns)
python3 - "$T/cr.tsv.gz" <<'PY'
import gzip, sys
cols = ["event_id", "gene_id", "te_id", "te_subfamily", "te_family", "te_class",
        "canonical", "chimera_type", "gene_te_distance", "n_samples", "total_reads"]
rows = [
    ["e_local", "GENE1", "TE_INTRON", "MTA", "ERVL-MaLR", "LTR", "yes", "te_initiated", "0", "2", "5"],
    ["e_far",   "GENE1", "TE_FAR",    "L1Md", "L1", "LINE", "no", "te_initiated", "350000", "1", "1"],
    ["e_trans", "GENE1", "TE_TRANS",  "L1Md", "L1", "LINE", "no", ".", "trans", "1", "1"],
    ["e_old",   "GENE1", "TE_OLD",    "B1", "Alu", "SINE", "no", "te_exonized", ".", "1", "1"],
]
with gzip.open(sys.argv[1], "wt") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
PY
if ! python3 workflow/scripts/chimera_evidence.py --junction "$T/cr.tsv.gz" \
      --cr-max-distance 200000 --out "$T/cand.tsv.gz" > "$T/ev.log" 2>&1; then
  echo "ERROR: chimera_evidence.py failed"; cat "$T/ev.log"; FAIL=1
fi

[ "$FAIL" = 0 ] || exit $FAIL

python3 - "$T" <<'PY' || FAIL=1
import csv, gzip, sys
T = sys.argv[1]
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

intr = {(r["chrom"], r["intron_start"], r["intron_end"]) for r in rows(f"{T}/introns.tsv.gz")}
check(intr == {("chr1", "1201", "2000"), ("chr1", "2201", "3000")},
      f"annotated introns must be 1-based inclusive like SJ.out.tab; got {sorted(intr)}")
last = open(f"{T}/last_exons.bed").read().split()
check(last[:3] == ["chr1", "3000", "3500"],
      f"last_exons.bed must hold TX1's 3'-most exon as BED (0-based start); got {last[:6]}")

sj = {r["intron_start"]: r for r in rows(f"{T}/sj.tsv")}
check(sj["1201"]["chimera_type"] == "annotated_splice"
      and sj["1201"]["gtf_annotated_intron"] == "yes",
      f"annotated intron with a splice site in a TE must be annotated_splice; got "
      f"{sj['1201']['chimera_type']!r}/{sj['1201'].get('gtf_annotated_intron')!r}")
check(sj["2701"]["chimera_type"] == "te_initiated"
      and sj["2701"]["gtf_annotated_intron"] == "no",
      f"novel TE->exon junction must stay te_initiated; got "
      f"{sj['2701']['chimera_type']!r}/{sj['2701'].get('gtf_annotated_intron')!r}")

asm = {r["transcript_id"]: r["chimera_type"] for r in rows(f"{T}/asm.tsv")}
check(asm.get("A_utr") == "annotated_terminal_exon_embedded_te",
      f"TE in the gene's annotated last exon must be annotated_terminal_exon_embedded_te, "
      f"got {asm.get('A_utr')!r}")
check(asm.get("A_novel") == "te_terminated",
      f"novel TE last exon must stay te_terminated, got {asm.get('A_novel')!r}")

cand = {(r["gene_id"], r["te_id"]): r for r in rows(f"{T}/cand.tsv.gz")}
check(("GENE1", "TE_INTRON") in cand, "local chimeric-read event must make a candidate")
check(("GENE1", "TE_FAR") not in cand,
      "event > max_gene_te_distance must NOT make a candidate")
check(("GENE1", "TE_TRANS") not in cand, "trans event must NOT make a candidate")
check(("GENE1", "TE_OLD") in cand,
      "event without a distance (older table) must be kept, not silently dropped")
check(cand.get(("GENE1", "TE_INTRON"), {}).get("cr_gene_te_distance") == "0",
      "cr_gene_te_distance must be carried for kept pairs")
log = open(f"{T}/ev.log").read()
check("2 chimeric-read event(s) skipped" in log,
      f"chimera_evidence.py must log how many far/trans events it skipped; log: {log!r}")
sys.exit(0 if ok else 1)
PY

exit $FAIL
