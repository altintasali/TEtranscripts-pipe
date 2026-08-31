#!/usr/bin/env bash
# Guard 25: chimera.assembly default (on) plans the full screen
#
# Run on its own:   .tests/guards/25_chimera_assembly_default_on_plans_the_full_screen.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# --- chimera.assembly.enabled defaults to TRUE as of 0.10.1, and
# config/test.yaml matches the shipped default on purpose so CI
# exercises the DAG users actually get. Every chimera_assembly
# rule, including its dedicated STAR pass and (test.yaml sets
# write_igv_bed: true) the IGV export, must be planned as-is.
if ! snakemake --configfile config/test.yaml -n --cores 1 >"$T/assembly_on.log" 2>&1; then
  echo "ERROR: default config must dry-run"; tail -40 "$T/assembly_on.log"; FAIL=1
else
  for rule in star_align_for_assembly stringtie_assemble stringtie_merge \
              stringtie_requantify chimera_assembly_classify chimera_assembly_quantify \
              chimera_assembly_cross_evidence chimera_assembly_summary_mqc chimera_assembly_igv_bed; do
    if ! grep -q "$rule" "$T/assembly_on.log"; then
      echo "ERROR: $rule not planned with the default chimera.assembly.enabled: true"
      tail -40 "$T/assembly_on.log"; FAIL=1
    fi
  done
fi

exit $FAIL
