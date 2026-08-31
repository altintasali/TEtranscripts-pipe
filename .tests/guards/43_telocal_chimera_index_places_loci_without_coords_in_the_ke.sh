#!/usr/bin/env bash
# Guard 43: TElocal chimera index places loci without coords in the key
#
# Run on its own:   .tests/guards/43_telocal_chimera_index_places_loci_without_coords_in_the_ke.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- A TElocal cntTable key is
# transcript_id:gene_id:family_id:class_id and carries NO
# coordinates unless the TE GTF's transcript_id happens to be a
# coordinate string (mghlab prebuilt annotations). The index used
# to parse coordinates out of the key, so on an ordinary rmsk-
# derived GTF (L1PA2_dupN) EVERY row was skipped, the index came
# out empty, and telocal_active was "no" for every junction --
# silently, with no error. Coordinates now come from
# telocal_locations.bed, and an unplaceable index is fatal.
mkdir -p "$T/tlx"
printf 'gene/TE\t/p/x.bam\nLx5b_dup2228:Lx5b:L1:LINE\t42\nMT2B1_dup5641:MT2B1:ERVL:LTR\t7\nENSMUSG00000000001\t900\n' \
  | gzip -c > "$T/tlx/s1.cntTable.gz"
printf 'chr2\t3000\t3500\tLx5b_dup2228:Lx5b:L1:LINE\t.\t+\nchr2\t9000\t9400\tMT2B1_dup5641:MT2B1:ERVL:LTR\t.\t-\n' \
  > "$T/tlx/loc.bed"
if ! python3 workflow/scripts/build_chimera_telocal_index.py \
    --telocal-tables "$T/tlx/s1.cntTable.gz" \
    --locations "$T/tlx/loc.bed" --out "$T/tlx/idx.pkl.gz" \
    > "$T/tlx/build.log" 2>&1; then
  echo "ERROR: index build failed"; cat "$T/tlx/build.log"; FAIL=1
fi
# ...and with no usable coordinates it must FAIL, not quietly
# produce an index that makes every locus look unexpressed.
printf 'chrX\t1\t2\tsomething_else\t.\t+\n' > "$T/tlx/empty.bed"
if python3 workflow/scripts/build_chimera_telocal_index.py \
    --telocal-tables "$T/tlx/s1.cntTable.gz" \
    --locations "$T/tlx/empty.bed" --out "$T/tlx/bad.pkl.gz" \
    > "$T/tlx/bad.log" 2>&1; then
  echo "ERROR: an unplaceable TElocal index must be fatal"; FAIL=1
elif ! grep -q "rows could be given genomic coordinates" "$T/tlx/bad.log"; then
  echo "ERROR: wrong message for an unplaceable index"; cat "$T/tlx/bad.log"; FAIL=1
fi
python3 - "$T/tlx" <<'PY' || FAIL=1
import sys
sys.path.insert(0, "workflow/scripts")
from chimera_telocal_index import TelocalIndex, classify_telocal, telocal_te_id
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
idx = TelocalIndex.load(f"{sys.argv[1]}/idx.pkl.gz")
hits = idx.overlapping("chr2", 3100, 3200)
check(len(hits) == 1, f"expected 1 overlapping locus, got {len(hits)}")
check(hits and hits[0][2] == "Lx5b_dup2228:Lx5b:L1:LINE", "wrong locus keyed")
check(hits and hits[0][5] == 42, "count did not survive the index")
check(hits and hits[0][3] == "L1", "family must resolve without coords in the key")
# both key shapes must classify, and gene rows must not
check(classify_telocal("Lx5b_dup2228:Lx5b:L1:LINE") == ("L1", "LINE"), "dupN key misclassified")
check(classify_telocal("chr1:5:9(L1PA2:+):L1PA2:L1:LINE") == ("L1", "LINE"), "coord key misclassified")
check(classify_telocal("ENSMUSG00000000001") == (None, None), "gene row must classify as None")
# te_id must be the INSERTION, which is what te.bed and the chimera tables use
check(telocal_te_id("Lx5b_dup2228:Lx5b:L1:LINE") == "Lx5b_dup2228", "te_id must be the insertion")
check(telocal_te_id("chr1:5:9(L1PA2:+):L1PA2:L1:LINE") == "chr1:5:9(L1PA2:+)", "coord-style te_id wrong")
sys.exit(0 if ok else 1)
PY

exit $FAIL
