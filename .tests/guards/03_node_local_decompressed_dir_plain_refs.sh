#!/usr/bin/env bash
# Guard 3: node-local decompressed_dir, plain refs
#
# Run on its own:   .tests/guards/03_node_local_decompressed_dir_plain_refs.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- node-local decompressed_dir is fine with plain refs ------------
# (no gunzip_reference jobs run, so the setting is dead config)
#
# telocal disabled here: config/test.yaml's telocal.locind is empty
# (auto-build), and telocal_locind now gzips its own output too (see
# rules/telocal.smk) -- decompressing THAT also goes through
# gunzip_reference, which would trip the node-local check below for a
# reason unrelated to what this guard tests (ref.fasta/gtf/te_gtf plain
# vs .gz). Guards 4/5/6 (fixture_gz_refs) already cover gunzip_reference
# + node-local decompressed_dir with telocal left at its default.
sed 's|  sjdb_overhang: 49|  sjdb_overhang: 49\n  decompressed_dir: /tmp/rnaseq-decompressed|' config/test.yaml > "$T/tmpdir_plain.yaml"
printf '\ntelocal:\n  enabled: false\n' >> "$T/tmpdir_plain.yaml"
if ! snakemake --configfile "$T/tmpdir_plain.yaml" -n --cores 1 >"$T/tmpdir_plain.log" 2>&1; then
  echo "ERROR: node-local decompressed_dir with plain refs must still dry-run"; FAIL=1
fi

exit $FAIL
