#!/usr/bin/env bash
# Guard 34: chimeric breakpoints resolve on BOTH strands at tolerance 0
#
# Run on its own:   .tests/guards/34_chimeric_breakpoints_resolve_on_both_strands_at_tolerance.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# STAR's Chimeric.out.junction reports the first base of the DONOR'S
# INTRON (col 2) and the last base of the ACCEPTOR'S INTRON (col 5)
# -- not the aligned segment edges. The aligned base is one step
# along the transcript, which flips sign with strand. The parser
# used to test an off-centre window around the raw breakpoint,
# which covered bp and bp+1 but never bp-1: acceptors landed inside
# it and donors fell just outside, so on real data donor
# breakpoints hit exons ~10x less often than acceptors and the
# gene_to_te / te_to_gene ratio was skewed.
#
# Both junctions below are genuine gene<->TE chimeras with
# breakpoints exactly on the feature edges. Before the fix the '+'
# case came back other_to_te (gene missed) and the '-' case
# gene_to_other (TE missed).
mkdir -p "$T/bpt"
printf 'chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";\n' > "$T/bpt/genes.gtf"
printf 'chr1\trmsk\texon\t3000\t3500\t.\t+\t.\tgene_id "L1"; transcript_id "L1_dup1"; family_id "L1"; class_id "LINE";\n' > "$T/bpt/te.gtf"
if ! python3 workflow/scripts/annotation_to_bed.py --gtf "$T/bpt/genes.gtf" \
      --te-gtf "$T/bpt/te.gtf" --outdir "$T/bpt/ref" > "$T/bpt/bed.log" 2>&1; then
  echo "ERROR: annotation_to_bed failed"; cat "$T/bpt/bed.log"; FAIL=1
else
  # + : donor ends at exon end 2000 -> STAR 2001 ; acceptor starts at TE 3000 -> STAR 2999
  printf 'chr1\t2001\t+\tchr1\t2999\t+\t1\t0\t0\tp1\t1\t50M\t1\t50M\n'  > "$T/bpt/j.junction"
  # - : mirrored. donor ends at 1000 -> STAR 999 ; acceptor starts at 3500 -> STAR 3501
  printf 'chr1\t999\t-\tchr1\t3501\t-\t1\t0\t0\tm1\t1\t50M\t1\t50M\n' >> "$T/bpt/j.junction"
  if ! python3 workflow/scripts/classify_chimera_reads.py \
        --junctions "$T/bpt/j.junction" --genes "$T/bpt/ref/genes.bed" \
        --exons "$T/bpt/ref/exons.bed" --te "$T/bpt/ref/te.bed" \
        --sample s --breakpoint-tolerance 0 --library-strandedness no \
        --out "$T/bpt/out.tsv.gz" > "$T/bpt/parse.log" 2>&1; then
    echo "ERROR: classify_chimera_reads.py failed"; cat "$T/bpt/parse.log"; FAIL=1
  else
    N=$(zcat "$T/bpt/out.tsv.gz" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["direction"]=="gene_to_te" && $h["gene_id"]=="G1" && $h["te_id"]=="L1_dup1"{n++} END{print n+0}')
    if [ "$N" -ne 2 ]; then
      echo "ERROR: expected 2 gene_to_te calls (one per strand), got $N"
      echo "  breakpoints sit exactly on the feature edges, so both must resolve"
      zcat "$T/bpt/out.tsv.gz" | cut -f5,8,15,17,19 | head; FAIL=1
    fi
  fi
fi

exit $FAIL
