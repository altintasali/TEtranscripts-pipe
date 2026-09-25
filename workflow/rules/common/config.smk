# Imports, config validation, per-rule resource lookup.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

import gzip
import os
import re
import tempfile

import pandas as pd
import yaml
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
# Renamed in 0.12.0: chimera.junction -> chimera.reads (itself since renamed
# again, see the next block). This runs BEFORE validate(), because the schema
# is additionalProperties: false and would otherwise reject the old key with
# a bare jsonschema error naming "junction" and nothing else -- which every
# existing user would hit, with no hint that the key moved. Fail fast, and
# say what to do about it.
if isinstance(config.get("chimera"), dict) and "junction" in config["chimera"]:
    raise WorkflowError(
        "config key 'chimera.junction' was renamed to 'chimera.reads' in "
        "0.12.0, and that key has itself since been renamed again -- see "
        "the 'chimera.reads' error below once you've made this rename.\n\n"
        "In your config file, rename:\n"
        "    chimera:\n"
        "      junction:      ->      reads:\n\n"
        "Everything nested under it is unchanged.\n\n"
        "The read screen's outputs also moved, from "
        "results/chimera/read_evidence/ to results/chimera/reads/ (and the "
        "assembly screen's from transcript_evidence/ to assembly/), so that "
        "stage re-runs once. Delete the old directories when you are happy "
        "with the new run."
    )

# Renamed: chimera.reads -> chimera.chimeric_reads. "Reads" was ambiguous --
# the SJ.out.tab screen is also read-derived -- and this screen's own report
# label ("Chimera (junction)"), config key ("reads") and data columns
# ("junction_*") were three different words for the same thing. Same
# fail-fast-before-validate() reasoning as the chimera.junction check above.
if isinstance(config.get("chimera"), dict) and "reads" in config["chimera"]:
    raise WorkflowError(
        "config key 'chimera.reads' was renamed to 'chimera.chimeric_reads'.\n"
        "In your config file, rename:\n"
        "    chimera:\n"
        "      reads:      ->      chimeric_reads:\n\n"
        "Everything nested under it is unchanged, except "
        "require_canonical_junction -> require_canonical (now matching "
        "chimera.splice_junctions.require_canonical's name).\n\n"
        "The screen's outputs also moved, from results/chimera/reads/ to "
        "results/chimera/chimeric_reads/, so that stage re-runs once. Delete "
        "the old directory when you are happy with the new run."
    )

# Renamed: chimera.sj_junctions -> chimera.splice_junctions. "sj_junctions"
# repeated itself ("SJ" already means splice junction), and the abbreviated
# "sj" was hard to read as an identifier elsewhere in the config/code: the
# spelled-out form is now used consistently everywhere except the already-
# compact report display text ("Chimera (SJ)", "SJ - ...").
if isinstance(config.get("chimera"), dict) and "sj_junctions" in config["chimera"]:
    raise WorkflowError(
        "config key 'chimera.sj_junctions' was renamed to "
        "'chimera.splice_junctions'.\n"
        "In your config file, rename:\n"
        "    chimera:\n"
        "      sj_junctions:      ->      splice_junctions:\n\n"
        "Everything nested under it is unchanged.\n\n"
        "The screen's outputs also moved, from results/chimera/sj/ to "
        "results/chimera/splice_junctions/, so that stage re-runs once. "
        "Delete the old directory when you are happy with the new run."
    )

validate(config, schema="../../schemas/config.schema.yaml")

V = config["versions"]

# -----------------------------------------------------------------------------
# Per-rule compute resources (threads/mem_mb/runtime), loaded from
# workflow/default-config/resources.yaml (built-in defaults) and optionally
# overridden by input/resources.yaml if you create one. A rule (or a
# missing key within its entry) not
# present there falls back to a small conservative default instead of
# failing -- see the "HPC / SLURM" section of the README for how these feed
# into cluster execution.
# -----------------------------------------------------------------------------
RESOURCES = config.get("resources", {})
_RESOURCE_DEFAULTS = {"threads": 1, "mem_mb": 4000, "runtime": 60}


def get_resources(rule_name):
    """Return {threads, mem_mb, runtime} for a rule name."""
    return {**_RESOURCE_DEFAULTS, **RESOURCES.get(rule_name, {})}


def get_scaled_mem_mb(rule_name):
    """Return *mem_mb* for *rule_name*, scaled by sample count.

    Rules that accumulate per-sample data in memory (chimera_chimeric_reads_counts,
    chimera_chimeric_reads_sample_qc_transform, tecount_counts, …) declare a ``mem_per_sample``
    key in ``resources.yaml`` on top of the base ``mem_mb``.  This helper
    computes ``base + per_sample × len(SAMPLES)`` so the SLURM allocation
    grows automatically with the experiment size.  Rules without
    ``mem_per_sample`` return the plain ``mem_mb`` value unchanged.
    """
    res = get_resources(rule_name)
    per_sample = RESOURCES.get(rule_name, {}).get("mem_per_sample", 0)
    return res["mem_mb"] + per_sample * len(SAMPLES)
