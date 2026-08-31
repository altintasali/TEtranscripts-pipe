#!/usr/bin/env bash
# Guard 44: index pickles are written atomically
#
# Run on its own:   .tests/guards/44_index_pickles_are_written_atomically.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- A real run died with a bare "EOFError: Ran out of input" from
# inside pickle.load: the chimera TElocal index on disk was 0 bytes,
# left by an earlier build that died mid-write. gzip and pickle have
# no length header, so nothing notices until a downstream job reads
# it, and the error names the reader rather than the failed write.
# Both index writers now go via a temp file + os.replace, so the
# final path is either absent or complete, and a corrupt one is
# reported with the file, its size, and the fix.
mkdir -p "$T/atomic"
printf 'gene/TE\t/p/x.bam\nLx5b_dup1:Lx5b:L1:LINE\t42\n' \
  | gzip -c > "$T/atomic/s.cntTable.gz"
printf 'chr2\t3000\t3500\tLx5b_dup1:Lx5b:L1:LINE\t.\t+\n' > "$T/atomic/loc.bed"
if ! python3 workflow/scripts/build_chimera_telocal_index.py \
    --telocal-tables "$T/atomic/s.cntTable.gz" \
    --locations "$T/atomic/loc.bed" --out "$T/atomic/idx.pkl.gz" \
    > "$T/atomic/build.log" 2>&1; then
  echo "ERROR: index build failed"; cat "$T/atomic/build.log"; FAIL=1
fi
python3 - "$T/atomic" <<'PY' || FAIL=1
import os, pickle, sys
sys.path.insert(0, "workflow/scripts")
from chimera_telocal_index import TelocalIndex
T = sys.argv[1]
ok = True
def check(cond, msg):
    global ok
    if not cond:
        print("ERROR:", msg); ok = False
idx = TelocalIndex.load(f"{T}/idx.pkl.gz")
check(sum(len(c["key"]) for c in idx.chroms.values()) == 1, "good index must still load")
# an empty file is exactly what the real failure had on disk
open(f"{T}/empty.pkl.gz", "wb").close()
try:
    TelocalIndex.load(f"{T}/empty.pkl.gz")
    check(False, "a 0-byte index must not load")
except SystemExit as e:
    msg = str(e)
    check("not a readable TElocal index" in msg, "unhelpful message for a corrupt index")
    check("0 bytes" in msg, "message must report the file size")
    check("rm " in msg, "message must say how to recover")
# a truncated gzip stream must be caught too, not just an empty file
data = open(f"{T}/idx.pkl.gz", "rb").read()
open(f"{T}/trunc.pkl.gz", "wb").write(data[: len(data) // 2])
try:
    TelocalIndex.load(f"{T}/trunc.pkl.gz")
    check(False, "a truncated index must not load")
except SystemExit:
    pass
# a write killed part-way must leave NO output and no temp file behind
real = pickle.dump
pickle.dump = lambda *a, **k: (_ for _ in ()).throw(KeyboardInterrupt("killed"))
try:
    idx.save(f"{T}/killed.pkl.gz")
    check(False, "interrupted save should propagate")
except BaseException:
    pass
finally:
    pickle.dump = real
check(not os.path.exists(f"{T}/killed.pkl.gz"), "interrupted save left a partial output")
check(not [f for f in os.listdir(T) if ".tmp." in f], "interrupted save left a temp file")
sys.exit(0 if ok else 1)
PY

exit $FAIL
