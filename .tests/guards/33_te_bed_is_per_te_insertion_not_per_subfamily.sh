#!/usr/bin/env bash
# Guard 33: te.bed is per TE INSERTION, not per subfamily
#
# Run on its own:   .tests/guards/33_te_bed_is_per_te_insertion_not_per_subfamily.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# annotation_to_bed.py used to key TE loci by gene_id -- the
# SUBFAMILY in the TEtranscripts GTF convention -- and take
# min(start)/max(end), collapsing every copy of a subfamily into
# one bounding box. On a real annotation that produced single
# intervals spanning ~192 Mb, so the breakpoint-overlap test could
# not fail and essentially every junction was called a TE hit.
# The bundled fixture has ONE TE, so only a multi-instance GTF
# catches this.
mkdir -p "$T/tebed"
printf 'chr1\trmsk\texon\t1000\t1200\t.\t+\t.\tgene_id "L1PA2"; transcript_id "L1PA2_dup1"; family_id "L1"; class_id "LINE";\n' > "$T/tebed/te.gtf"
printf 'chr1\trmsk\texon\t900000\t900200\t.\t+\t.\tgene_id "L1PA2"; transcript_id "L1PA2_dup2"; family_id "L1"; class_id "LINE";\n' >> "$T/tebed/te.gtf"
# two rows sharing a transcript_id = one fragmented insertion -> must MERGE
printf 'chr1\trmsk\texon\t7000\t7100\t.\t+\t.\tgene_id "L1PA2"; transcript_id "L1PA2_dup3"; family_id "L1"; class_id "LINE";\n' >> "$T/tebed/te.gtf"
printf 'chr1\trmsk\texon\t7150\t7300\t.\t+\t.\tgene_id "L1PA2"; transcript_id "L1PA2_dup3"; family_id "L1"; class_id "LINE";\n' >> "$T/tebed/te.gtf"
printf 'chr1\tsrc\texon\t2000\t2500\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";\n' > "$T/tebed/genes.gtf"
if ! python3 workflow/scripts/annotation_to_bed.py --gtf "$T/tebed/genes.gtf" \
      --te-gtf "$T/tebed/te.gtf" --outdir "$T/tebed/out" > "$T/tebed.log" 2>&1; then
  echo "ERROR: annotation_to_bed.py failed"; cat "$T/tebed.log"; FAIL=1
else
  N=$(wc -l < "$T/tebed/out/te.bed")
  if [ "$N" -ne 3 ]; then
    echo "ERROR: te.bed should have 3 rows (dup1, dup2, merged dup3), got $N"
    cat "$T/tebed/out/te.bed"; FAIL=1
  fi
  # the regression itself: no interval may span the two distant copies
  if awk -F'\t' '$3 - $2 > 100000 {found=1} END{exit !found}' "$T/tebed/out/te.bed"; then
    echo "ERROR: te.bed has a >100kb interval -- subfamily bounding-box bug is back"
    cat "$T/tebed/out/te.bed"; FAIL=1
  fi
  # fragmented insertion merged into one 7000-7300 interval
  if ! awk -F'\t' '$4=="L1PA2_dup3" && $2==6999 && $3==7300 {found=1} END{exit !found}' "$T/tebed/out/te.bed"; then
    echo "ERROR: fragments sharing a transcript_id were not merged"
    cat "$T/tebed/out/te.bed"; FAIL=1
  fi
  # subfamily kept in col 9 (TEcount cross-reference depends on it)
  if ! awk -F'\t' '$4=="L1PA2_dup1" && $9=="L1PA2" {found=1} END{exit !found}' "$T/tebed/out/te.bed"; then
    echo "ERROR: te.bed col 9 should carry the subfamily"; cat "$T/tebed/out/te.bed"; FAIL=1
  fi
fi

exit $FAIL
