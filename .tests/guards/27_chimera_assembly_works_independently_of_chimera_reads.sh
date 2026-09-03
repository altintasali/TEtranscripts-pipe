#!/usr/bin/env bash
# Guard 27: chimera.assembly works independently of chimera.reads
#
# Run on its own:   .tests/guards/27_chimera_assembly_works_independently_of_chimera_reads.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- assembly on, reads screen OFF: annotation_to_bed (moved to
# ref.smk specifically so this works) must still run, and
# chimera_assembly_cross_evidence (which needs the reads
# screen's output) must NOT be planned.
# (derives from config/test.yaml, where assembly is already on by
# default -- there is no separate assembly_on.yaml to build from.)
sed '/^  reads:$/{n;s/enabled: true/enabled: false/}' \
    config/test.yaml > "$T/assembly_only.yaml"
if ! grep -qE "^    enabled: false" "$T/assembly_only.yaml"; then
  echo "ERROR: guard 27 sed did not turn the reads screen off"
  sed -n '/^  reads:/,+2p' "$T/assembly_only.yaml"; FAIL=1
elif ! snakemake --configfile "$T/assembly_only.yaml" -n --cores 1 >"$T/assembly_only.log" 2>&1; then
  echo "ERROR: assembly-only config must dry-run"; tail -40 "$T/assembly_only.log"; FAIL=1
elif ! grep -q "annotation_to_bed" "$T/assembly_only.log"; then
  echo "ERROR: annotation_to_bed not planned with chimera.reads disabled "
  echo "(chimera_assembly_classify needs it -- see ref.smk)"
  tail -40 "$T/assembly_only.log"; FAIL=1
elif grep -q "chimera_assembly_cross_evidence" "$T/assembly_only.log"; then
  echo "ERROR: chimera_assembly_cross_evidence planned with chimera.reads disabled"
  tail -40 "$T/assembly_only.log"; FAIL=1
elif ! grep -q "chimera_assembly_classify" "$T/assembly_only.log"; then
  echo "ERROR: chimera_assembly_classify not planned with chimera.reads disabled"
  tail -40 "$T/assembly_only.log"; FAIL=1
fi

exit $FAIL
