#!/usr/bin/env bash
# Guard 53: the standalone candidates explorer shows every row and no score
#
# Run on its own:   .tests/guards/53_candidates_explorer_shows_every_row_no_score.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- chimera_candidates_explorer.R is the "every one of tens of thousands of
# rows" companion to chimera_candidates_table_mqc.py's top-N MultiQC table
# (embedding the full catalogue there runs the report from 2.3 MB to 83 MB
# and MultiQC's own build time 10x, measured directly). It also adds two
# columns candidates.tsv.gz does not carry: a ready-to-paste IGV locus per
# gene and per TE, joined from genes.bed/te.bed (both keyed on the same
# gene_id/te_id candidates.tsv.gz uses) -- BED is 0-based half-open, IGV's
# locus box is 1-based inclusive, so this pins that off-by-one exactly.
#
# Same non-ranking stance as guards 36/50: every evidence column is shown
# as-is and no combined/weighted score is ever computed here.
if ! command -v Rscript >/dev/null 2>&1 || ! Rscript -e 'library(DT); library(htmlwidgets)' >/dev/null 2>&1; then
  echo "[guard 53] SKIP: R/DT/htmlwidgets not on PATH in this environment " \
       "(they live in the generated candidates_explorer env, only present " \
       "under --sdm conda or after installing workflow/environment.yaml's " \
       "r-dt/r-htmlwidgets/pandoc directly)"
  exit 0
fi

mkdir -p "$T/exp"
cols="gene_id\tte_id\tte_subfamily\tte_family\tte_class\tfound_by\tevidence\tn_evidence\tjunction_events\tjunction_reads\tjunction_max_samples\tjunction_canonical\tjunction_chimera_types\ttelocal_active\ttelocal_count\tassembly_transcripts\tassembly_chimera_types\tassembly_strand_match\tassembly_transcript_ids"
{
  printf "%b\n" "$cols"
  # resolved-symbol row (Gapdh), telocal ran and found reads
  printf 'ENSMUSG00000057666\tL1PA2_dup1\tL1PA2\tL1\tLINE\tboth\tcanonical,multi_sample,both_screens,assembly_strand_match,telocal_expressed\t5\t2\t60\t3\tyes\tte_terminated\tyes\t42\t2\tte_terminated\tyes\tMSTRG.1.1\n'
  # gene_id-fallback row (no symbol), telocal NEVER RAN ("." must stay blank, not 0)
  printf 'ENSMUSG99999999999\tAluY_dup9\tAluY\tAlu\tSINE\tassembly\t.\t0\t0\t0\t0\tno\t.\t.\t.\t1\tte_exonized\t.\tMSTRG.2.1\n'
} | gzip -c > "$T/exp/candidates.tsv.gz"

{ printf 'gene_id\tgene_name\n'; printf 'ENSMUSG00000057666\tGapdh\n'; } \
  | gzip -c > "$T/exp/gene_id_to_name.tsv.gz"

# BED: chrom start end id score strand [family class subfamily] -- 0-based
# start, no header (annotation_to_bed.py's convention).
printf 'chr1\t999\t2000\tENSMUSG00000057666\t.\t+\nchr2\t4999\t6000\tENSMUSG99999999999\t.\t-\n' \
  > "$T/exp/genes.bed"
printf 'chr1\t2099\t2400\tL1PA2_dup1\t.\t+\tL1\tLINE\tL1PA2\nchr2\t7099\t7300\tAluY_dup9\t.\t-\tAlu\tSINE\tAluY\n' \
  > "$T/exp/te.bed"

if ! Rscript workflow/scripts/chimera_candidates_explorer.R \
      "$T/exp/candidates.tsv.gz" "$T/exp/gene_id_to_name.tsv.gz" \
      "$T/exp/genes.bed" "$T/exp/te.bed" "$T/exp/out.html" \
      > "$T/exp/log" 2>&1; then
  echo "ERROR: chimera_candidates_explorer.R failed"; cat "$T/exp/log"; FAIL=1
else
  size=$(wc -c < "$T/exp/out.html" | tr -d ' ')
  if [ "$size" -lt 100000 ]; then
    echo "ERROR: output HTML is only $size bytes -- too small to be a real, self-contained DT bundle (expected > 100000)"
    FAIL=1
  fi
  if ! grep -qi "datatables" "$T/exp/out.html"; then
    echo "ERROR: no recognizable DataTables markup in the output"; FAIL=1
  fi
  for h in "Gene" "TE insertion" "Gene locus" "TE locus" "TE subfamily" \
           "TE family" "TE class" "Found by" "Evidence flags" \
           "Evidence types" "Splice motif" "Samples" "Junction events" \
           "Reads" "Junction chimera types" "TElocal active" \
           "TE locus reads" "Assembly transcripts" "Assembly chimera types" \
           "Strand match" "Assembly transcript IDs"; do
    if ! grep -qF "\"$h\"" "$T/exp/out.html"; then
      echo "ERROR: expected column header missing: $h"; FAIL=1
    fi
  done
  # BED "chr1 999 2000" -> IGV locus "chr1:1000-2000": pins the +1 exactly.
  if ! grep -qF "chr1:1000-2000" "$T/exp/out.html"; then
    echo "ERROR: gene locus off-by-one wrong -- expected chr1:1000-2000 (BED start 999 + 1)"
    FAIL=1
  fi
  if ! grep -qF "chr1:2100-2400" "$T/exp/out.html"; then
    echo "ERROR: TE locus off-by-one wrong -- expected chr1:2100-2400 (BED start 2099 + 1)"
    FAIL=1
  fi
  # both symbol-resolution paths must survive to the rendered page.
  if ! grep -qF "Gapdh" "$T/exp/out.html"; then
    echo "ERROR: resolved gene symbol (Gapdh) missing"; FAIL=1
  fi
  if ! grep -qF "ENSMUSG99999999999" "$T/exp/out.html"; then
    echo "ERROR: gene_id fallback (no symbol available) missing"; FAIL=1
  fi
  # anti-regression, same spirit as guards 36/50: no combined/weighted score
  # column may ever appear here, only the raw per-column values.
  if grep -oE '"[^"]*[Cc]onfidence[^"]*"|"[Rr]ank"|"[Ss]core"' "$T/exp/out.html" \
       | grep -v '"score"'; then
    echo "ERROR: a Confidence/Rank/Score-like column header was found -- this table must never carry a combined ranking"
    FAIL=1
  fi
fi

exit $FAIL
