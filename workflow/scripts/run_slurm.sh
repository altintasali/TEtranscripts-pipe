#!/usr/bin/env bash
# =============================================================================
# Convenience launcher for running the workflow on a SLURM cluster.
#
# Usage:
#   ./workflow/scripts/run_slurm.sh [snakemake options...]
#
# Equivalent to:
#   snakemake --workflow-profile workflow/profiles/slurm <your options>
#
# All extra arguments are passed straight through to snakemake, e.g.:
#   ./workflow/scripts/run_slurm.sh --configfile config/test.yaml  # smoke test
#   ./workflow/scripts/run_slurm.sh --cores 64                     # limit CPU cores
#   ./workflow/scripts/run_slurm.sh -n                             # dry-run
#
# The SLURM profile (workflow/profiles/slurm/config.yaml) submits every job
# with sbatch and sizes each job from config/resources.yaml; see the README
# "HPC / SLURM" section for overriding resources and the account/partition.
# =============================================================================
set -euo pipefail

# Always run from the repo root, wherever this script is invoked from, so
# that the relative paths in the Snakefile (config/, input/, results/, logs/)
# and in the SLURM profile resolve the same way as a local run.
cd "$(dirname "$0")/../.."

exec snakemake --workflow-profile workflow/profiles/slurm "$@"
