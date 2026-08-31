#!/usr/bin/env bash
# Guard 31: telocal.indexer: legacy opts into TElocal's own build()
#
# Run on its own:   .tests/guards/31_telocal_indexer_legacy_opts_into_telocal_s_own_build.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

sed '/^  locind: ""$/a\  indexer: legacy' config/test.yaml > "$T/indexer_legacy.yaml"
if ! snakemake --configfile "$T/indexer_legacy.yaml" -n -p --cores 1 --until telocal_locind >"$T/indexer_legacy.log" 2>&1; then
  echo "ERROR: telocal.indexer: legacy must dry-run"; tail -40 "$T/indexer_legacy.log"; FAIL=1
elif ! grep -q -- "--indexer legacy" "$T/indexer_legacy.log"; then
  echo "ERROR: telocal.indexer: legacy did not produce --indexer legacy"
  tail -40 "$T/indexer_legacy.log"; FAIL=1
fi

# --- both indexers must produce equivalent output on the real
# small test fixture (metadata equality; the exhaustive spatial-
# query verification lives in build_telocal_index.py's own
# commit history, run manually against real TElocal_Toolkit).
if ! python3 workflow/scripts/build_telocal_index.py \
      --gtf .tests/resources/te_annotation.gtf --indexer fast \
      --out "$T/idx_fast.locInd" > "$T/idx_fast.log" 2>&1; then
  echo "ERROR: build_telocal_index.py --indexer fast failed"; cat "$T/idx_fast.log"; FAIL=1
fi
if ! python3 workflow/scripts/build_telocal_index.py \
      --gtf .tests/resources/te_annotation.gtf --indexer legacy \
      --out "$T/idx_legacy.locInd" > "$T/idx_legacy.log" 2>&1; then
  echo "ERROR: build_telocal_index.py --indexer legacy failed"; cat "$T/idx_legacy.log"; FAIL=1
fi
if [ -f "$T/idx_fast.locInd" ] && [ -f "$T/idx_legacy.locInd" ]; then
  if ! python3 -c "import pickle,sys; a=pickle.load(open('$T/idx_fast.locInd','rb')); b=pickle.load(open('$T/idx_legacy.locInd','rb')); sys.exit(0 if (a.getNames()==b.getNames() and a.getElements()==b.getElements() and a._length==b._length) else 1)"; then
    echo "ERROR: fast and legacy indexers produced different output"; FAIL=1
  else
    echo "fast/legacy indexer output: IDENTICAL"
  fi
fi

exit $FAIL
