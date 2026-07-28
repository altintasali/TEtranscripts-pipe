import itertools
import os

import pandas as pd
import yaml
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

V = config["versions"]

# -----------------------------------------------------------------------------
# Load & validate sample sheet
# columns: sample, fastq_1, fastq_2 (optional), condition (optional)
# -----------------------------------------------------------------------------
samples = (
    pd.read_csv(config["samples"], dtype=str, comment="#")
    .set_index("sample", drop=False)
    .sort_index()
)
# jsonschema's "null" type only matches Python None, not pandas/NumPy NaN, so
# empty optional cells (fastq_2, condition) need to be converted explicitly.
samples = samples.where(pd.notnull(samples), None)
validate(samples, schema="../schemas/samples.schema.yaml")

SAMPLES = list(samples["sample"])
HAS_CONDITION = "condition" in samples.columns and samples["condition"].notna().all()


def _is_paired(sample):
    """A sample is paired-end if it has a non-empty fastq_2 value.

    This is resolved per-sample (not globally), so a single sample sheet can
    freely mix paired-end and single-end samples.
    """
    if "fastq_2" not in samples.columns:
        return False
    val = samples.loc[sample, "fastq_2"]
    return pd.notna(val) and str(val).strip() != ""


def star_input(wildcards):
    """Return {fq1: ..., [fq2: ...]} paths for a sample."""
    d = {"fq1": samples.loc[wildcards.sample, "fastq_1"]}
    if _is_paired(wildcards.sample):
        d["fq2"] = samples.loc[wildcards.sample, "fastq_2"]
    return d


def star_reads_param(wildcards):
    """Space-separated --readFilesIn argument (fq1 [fq2])."""
    d = star_input(wildcards)
    reads = [d["fq1"]]
    if "fq2" in d:
        reads.append(d["fq2"])
    return " ".join(reads)


def star_read_command_param(wildcards):
    """--readFilesCommand zcat if the fastq files are gzip-compressed."""
    d = star_input(wildcards)
    if str(d["fq1"]).endswith(".gz"):
        return "--readFilesCommand zcat"
    return ""


def strandedness_input(wildcards):
    """Dependency on the per-sample auto-detected strandedness call.

    Returns [] (no dependency) when strandedness.mode is fixed in the config,
    since no RSeQC/auto-detection step is run in that case.
    """
    if config["strandedness"]["mode"] == "auto":
        return f"results/rseqc/{wildcards.sample}/strandedness.txt"
    return []


def get_strandedness_param(wildcards, input):
    """Resolve the --stranded value (no/forward/reverse) for a sample.

    If strandedness.mode is "auto", read the value that
    workflow/scripts/determine_strandedness.py determined from the RSeQC
    infer_experiment.py output (the file is guaranteed to already exist
    because it is declared as a rule input alongside this params function).
    Otherwise, the fixed value from the config is used for every sample.
    """
    mode = config["strandedness"]["mode"]
    if mode == "auto":
        with open(input.strandedness) as fh:
            return fh.read().strip()
    return mode


# -----------------------------------------------------------------------------
# Pairwise contrasts, derived automatically from the "condition" column.
# If the sample sheet has no "condition" column, CONTRASTS is empty and the
# TEtranscripts differential-analysis step is skipped entirely.
# -----------------------------------------------------------------------------
def _build_contrasts():
    if not HAS_CONDITION:
        return {}
    conditions = sorted(samples["condition"].unique())
    contrasts = {}
    for control, treatment in itertools.combinations(conditions, 2):
        name = f"{treatment}_vs_{control}"
        contrasts[name] = {
            "treatment": list(samples.loc[samples["condition"] == treatment, "sample"]),
            "control": list(samples.loc[samples["condition"] == control, "sample"]),
        }
    return contrasts


CONTRASTS = _build_contrasts()


def contrast_strandedness_input(wildcards):
    """Dependency on the strandedness call of every sample in a contrast."""
    if config["strandedness"]["mode"] != "auto":
        return []
    c = CONTRASTS[wildcards.contrast]
    all_samples = list(c["treatment"]) + list(c["control"])
    return expand("results/rseqc/{sample}/strandedness.txt", sample=all_samples)


def get_contrast_strandedness_param(wildcards, input):
    """Resolve a single --stranded value shared by every sample in a contrast.

    TEtranscripts runs one DESeq2 analysis across all treatment/control BAMs
    at once, so it needs one strandedness value. If auto-detection disagrees
    between samples in the same contrast, this fails loudly rather than
    silently picking one -- that mismatch usually means samples were
    prepared with different library kits and shouldn't be pooled blindly.
    """
    mode = config["strandedness"]["mode"]
    if mode != "auto":
        return mode
    values = set()
    for f in input.strandedness:
        with open(f) as fh:
            values.add(fh.read().strip())
    if len(values) > 1:
        raise ValueError(
            f"Samples in contrast '{wildcards.contrast}' have inconsistent "
            f"auto-detected strandedness ({values}). Verify these samples "
            "were prepared with the same library protocol, or set "
            "strandedness.mode explicitly in config.yaml instead of 'auto'."
        )
    return values.pop()


def all_tecount_tables():
    return expand("results/tecount/{sample}/{sample}.cntTable", sample=SAMPLES)


def all_diffexp_outputs():
    return expand(
        "results/tetranscripts/{contrast}/{contrast}_sigdiff_gene_TE.txt",
        contrast=list(CONTRASTS.keys()),
    )


# -----------------------------------------------------------------------------
# Conda environments for every tool in the workflow, generated from
# config["versions"]. Editing a version string there and re-running is all
# that's needed -- Snakemake/conda resolves and downloads the matching build
# automatically (via `--sdm conda` / `--use-conda`), no manual install
# required. This intentionally does NOT rely on snakemake-wrapper-bundled
# environments, because a wrapper's pinned tag moves all of its tools'
# versions together -- it can't give you e.g. STAR 2.7.11b + samtools 1.21
# independently if that particular wrapper tag happens to pin a different
# samtools. Generating one small env per tool from explicit version strings
# is the only way to match an arbitrary external reference (e.g. a specific
# nf-core/rnaseq release's tool versions) exactly.
# -----------------------------------------------------------------------------
GENERATED_ENV_DIR = "workflow/envs/generated"
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)


def _write_env(name, dependencies):
    path = f"{GENERATED_ENV_DIR}/{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(
            {"channels": ["bioconda", "conda-forge"], "dependencies": dependencies},
            fh,
            sort_keys=False,
        )
    return path


STAR_ENV = _write_env("star", [f"star={V['star']}"])
SAMTOOLS_ENV = _write_env("samtools", [f"samtools={V['samtools']}"])
RSEQC_ENV = _write_env("rseqc", [f"rseqc={V['rseqc']}"])
MULTIQC_ENV = _write_env("multiqc", [f"multiqc={V['multiqc']}"])

TETRANSCRIPTS_ENV = _write_env(
    "tetranscripts",
    [
        f"tetranscripts={V['tetranscripts']}",
        f"bioconductor-deseq2={V['deseq2']}",
        "r-base",
    ],
)

UCSC_TOOLS_ENV = _write_env(
    "ucsc_tools",
    [
        f"ucsc-gtftogenepred={V['ucsc_gtftogenepred']}",
        f"ucsc-genepredtobed={V['ucsc_genepredtobed']}",
    ],
)
