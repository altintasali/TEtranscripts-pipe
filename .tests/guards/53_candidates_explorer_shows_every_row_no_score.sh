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
# locus box is 1-based inclusive, so this pins that off-by-one exactly. It
# also joins a real cohort-total read count per candidate from
# chimera_candidates_matrix_totals.py's pre-summed output, when the
# corresponding screen ran.
#
# Same non-ranking stance as guards 36/50: every evidence column is shown
# as-is and no combined/weighted score is ever computed here.
mkdir -p "$T/exp"
cols="gene_id\tte_id\tte_subfamily\tte_family\tte_class\tfound_by\tevidence\tn_evidence\tjunction_events\tjunction_reads\tjunction_max_samples\tjunction_canonical\tjunction_chimera_types\ttelocal_active\ttelocal_count\ttelocal_locus\tassembly_transcripts\tassembly_chimera_types\tassembly_strand_match\tassembly_transcript_ids"
{
  printf "%b\n" "$cols"
  # resolved-symbol row (Gapdh), telocal ran and found reads, both cohort
  # totals available (telocal_locus / assembly_transcript_ids match the
  # totals fixtures below)
  printf 'ENSMUSG00000057666\tL1PA2_dup1\tL1PA2\tL1\tLINE\tboth\tcanonical,multi_sample,both_screens,assembly_strand_match,telocal_expressed\t5\t2\t60\t3\tyes\tte_terminated\tyes\t42\tL1PA2_dup1:L1PA2:L1PA2fam:LINE\t2\tte_terminated\tyes\tMSTRG.1.1\n'
  # gene_id-fallback row (no symbol), telocal NEVER RAN ("." must stay blank,
  # not 0) and its assembly transcript has no matching totals row (blank too)
  printf 'ENSMUSG99999999999\tAluY_dup9\tAluY\tAlu\tSINE\tassembly\t.\t0\t0\t0\t0\tno\t.\t.\t.\t.\t1\tte_exonized\t.\tMSTRG.2.1\n'
} | gzip -c > "$T/exp/candidates.tsv.gz"

{ printf 'gene_id\tgene_name\n'; printf 'ENSMUSG00000057666\tGapdh\n'; } \
  | gzip -c > "$T/exp/gene_id_to_name.tsv.gz"

# BED: chrom start end id score strand [family class subfamily] -- 0-based
# start, no header (annotation_to_bed.py's convention).
printf 'chr1\t999\t2000\tENSMUSG00000057666\t.\t+\nchr2\t4999\t6000\tENSMUSG99999999999\t.\t-\n' \
  > "$T/exp/genes.bed"
printf 'chr1\t2099\t2400\tL1PA2_dup1\t.\t+\tL1\tLINE\tL1PA2\nchr2\t7099\t7300\tAluY_dup9\t.\t-\tAlu\tSINE\tAluY\n' \
  > "$T/exp/te.bed"

# chimera_candidates_matrix_totals.py's own output shape (key <TAB> total).
# Only the FIRST candidate row's keys are present, so the second row's
# cohort-total columns must resolve to blank (not 0 or an error).
printf 'key\ttotal\nL1PA2_dup1:L1PA2:L1PA2fam:LINE\t123.456\n' \
  > "$T/exp/telocal_totals.tsv"
printf 'key\ttotal\nMSTRG.1.1\t77.000\n' \
  > "$T/exp/assembly_totals.tsv"

# --- The rule's own shell step (not the R script) awk-filters genes.bed/
# te.bed down to just the candidates' gene_id/te_id before R ever sees them.
# genes.bed/te.bed are GENOME-WIDE (every gene/TE in the annotation, not just
# candidates) -- te.bed especially, millions of rows for a real mouse/human
# TE annotation. Loading the unfiltered file OOM-killed a real run (2.3 GB
# against a 2.4 GB request) even though a candidate-count-sized dev fixture
# never caught it. This exercises the EXACT awk command from
# workflow/rules/chimera_reads.smk (kept in sync by eye -- if that command
# changes, update this copy too) against a small stand-in for a genome-wide
# BED that mixes candidate and non-candidate ids, and pins that only the
# candidate rows survive. Needs no R, so it runs even where Rscript/DT are
# unavailable.
gzip -dc "$T/exp/candidates.tsv.gz" | tail -n +2 | cut -f1 | sort -u > "$T/exp/gene_ids.txt"
gzip -dc "$T/exp/candidates.tsv.gz" | tail -n +2 | cut -f2 | sort -u > "$T/exp/te_ids.txt"
{
  cat "$T/exp/genes.bed"
  # non-candidate genes that must be FILTERED OUT
  printf 'chr3\t0\t100\tENSMUSG_NOT_A_CANDIDATE_1\t.\t+\n'
  printf 'chr4\t0\t100\tENSMUSG_NOT_A_CANDIDATE_2\t.\t+\n'
} > "$T/exp/genome_wide_genes.bed"
{
  cat "$T/exp/te.bed"
  printf 'chr3\t0\t100\tL1_NOT_A_CANDIDATE\t.\t+\tL1\tLINE\tL1x\n'
} > "$T/exp/genome_wide_te.bed"
awk -F'\t' 'NR==FNR{ids[$1]=1; next} ($4 in ids)' \
  "$T/exp/gene_ids.txt" "$T/exp/genome_wide_genes.bed" > "$T/exp/filtered_genes.bed"
awk -F'\t' 'NR==FNR{ids[$1]=1; next} ($4 in ids)' \
  "$T/exp/te_ids.txt" "$T/exp/genome_wide_te.bed" > "$T/exp/filtered_te.bed"
n_genes=$(wc -l < "$T/exp/filtered_genes.bed" | tr -d ' ')
n_te=$(wc -l < "$T/exp/filtered_te.bed" | tr -d ' ')
if [ "$n_genes" != "2" ]; then
  echo "ERROR: gene BED filter kept $n_genes rows, expected exactly the 2 candidate genes -- non-candidate genome rows leaked through, or a candidate was dropped"
  FAIL=1
fi
if [ "$n_te" != "2" ]; then
  echo "ERROR: TE BED filter kept $n_te rows, expected exactly the 2 candidate TEs -- non-candidate genome rows leaked through, or a candidate was dropped"
  FAIL=1
fi
if grep -q "NOT_A_CANDIDATE" "$T/exp/filtered_genes.bed" "$T/exp/filtered_te.bed"; then
  echo "ERROR: a non-candidate genome-wide row survived filtering -- this is exactly what OOM-killed a real run"
  FAIL=1
fi

# --- chimera_candidates_matrix_totals.py itself: streams a matrix once and
# sums only the requested keys. Needs no R, runs even where Rscript/DT are
# unavailable.
printf 'transcript_id\ts1\ts2\nMSTRG.1.1\t30\t47\nMSTRG.9.9\t999\t999\n' \
  | gzip -c > "$T/exp/assembly_counts_matrix.tsv.gz"
printf 'MSTRG.1.1\n' > "$T/exp/assembly_keys.txt"
if ! python3 workflow/scripts/chimera_candidates_matrix_totals.py \
      --matrix "$T/exp/assembly_counts_matrix.tsv.gz" \
      --keys "$T/exp/assembly_keys.txt" --out "$T/exp/totals_check.tsv" \
      > "$T/exp/totals_check.log" 2>&1; then
  echo "ERROR: chimera_candidates_matrix_totals.py failed"; cat "$T/exp/totals_check.log"; FAIL=1
else
  if ! grep -qF "MSTRG.1.1	77.000" "$T/exp/totals_check.tsv"; then
    echo "ERROR: expected MSTRG.1.1 total 77.000 (30+47), got:"; cat "$T/exp/totals_check.tsv"; FAIL=1
  fi
  if grep -q "MSTRG.9.9" "$T/exp/totals_check.tsv"; then
    echo "ERROR: a key NOT in --keys leaked into the totals output -- this is exactly the genome-scale-payload risk the streaming filter exists to avoid"
    FAIL=1
  fi
fi

# --- ANTI-REGRESSION: the rule's shell runs under `set -euo pipefail`
# (conda-env activation implies it), and `grep -v '^\.$' | sort -u > file`
# exits 1 when EVERY row is "." -- pipefail picks up grep's non-zero exit
# even though sort (the actual last command) succeeds, aborting the whole
# `&&` chain. This is exactly what broke a real end-to-end run: a
# candidates.tsv.gz where every assembly_transcript_ids was "." (no
# assembly-screen candidates at all) killed chimera_candidates_explorer
# outright. Pins that the rule's own key-extraction lines (kept in sync by
# eye with chimera_reads.smk's _candidates_explorer_shell(), same caveat as
# the awk commands above) tolerate an all-"." column.
printf '.\n.\n.\n' | gzip -c > "$T/exp/all_dot_candidates.tsv.gz"
if ! bash -c '
  set -euo pipefail
  gzip -dc "$1" | grep -v "^\.\$" | sort -u > "$2" || true
' _ "$T/exp/all_dot_candidates.tsv.gz" "$T/exp/all_dot_keys.txt"; then
  echo "ERROR: the rule's key-extraction line aborts under pipefail when every row is '.' -- this is the exact bug that broke a real run"
  FAIL=1
elif [ -s "$T/exp/all_dot_keys.txt" ]; then
  echo "ERROR: expected an empty keys file when every candidate row is '.', got:"
  cat "$T/exp/all_dot_keys.txt"; FAIL=1
fi

if ! command -v Rscript >/dev/null 2>&1 || ! Rscript -e 'library(DT); library(htmlwidgets)' >/dev/null 2>&1; then
  echo "[guard 53] SKIP (R portion only): R/DT/htmlwidgets not on PATH in " \
       "this environment (they live in the generated candidates_explorer " \
       "env, only present under --sdm conda or after installing " \
       "workflow/environment.yaml's r-dt/r-htmlwidgets/pandoc directly). " \
       "The non-R checks above still ran."
  exit $FAIL
fi

if ! Rscript workflow/scripts/chimera_candidates_explorer.R \
      "$T/exp/candidates.tsv.gz" "$T/exp/gene_id_to_name.tsv.gz" \
      "$T/exp/genes.bed" "$T/exp/te.bed" \
      "$T/exp/telocal_totals.tsv" "$T/exp/assembly_totals.tsv" \
      "$T/exp/out.html" \
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
           "Evidence count" "Splice motif" "Chimeric junction samples" \
           "Junction events" "Junction reads (cohort total)" "TE type (reads)" \
           "TElocal active" \
           "TElocal reads (cohort total)" "Assembly transcript count" \
           "Assembly reads (cohort total)" "TE type (assembly)" \
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
  # the cohort-total joins: row 1 gets real numbers, row 2 (no matching
  # totals key) must stay blank, not silently become 0.
  if ! grep -qF "123.456" "$T/exp/out.html"; then
    echo "ERROR: TElocal cohort-total join missing (expected 123.456)"; FAIL=1
  fi
  if ! grep -qF "77" "$T/exp/out.html"; then
    echo "ERROR: Assembly cohort-total join missing (expected 77)"; FAIL=1
  fi
  # export buttons and the injected font, both additive UI requests.
  if ! grep -qi "buttons" "$T/exp/out.html"; then
    echo "ERROR: no Buttons extension (copy/csv/excel) markup found"; FAIL=1
  fi
  if ! grep -qF "Helvetica" "$T/exp/out.html"; then
    echo "ERROR: injected font-family CSS missing"; FAIL=1
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
