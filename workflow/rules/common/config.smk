# Imports, config validation, per-rule resource lookup.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

import gzip
import itertools
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

    Rules that accumulate per-sample data in memory (chimera_counts,
    sample_qc_transform, tecount_counts, …) declare a ``mem_per_sample``
    key in ``resources.yaml`` on top of the base ``mem_mb``.  This helper
    computes ``base + per_sample × len(SAMPLES)`` so the SLURM allocation
    grows automatically with the experiment size.  Rules without
    ``mem_per_sample`` return the plain ``mem_mb`` value unchanged.
    """
    res = get_resources(rule_name)
    per_sample = RESOURCES.get(rule_name, {}).get("mem_per_sample", 0)
    return res["mem_mb"] + per_sample * len(SAMPLES)
