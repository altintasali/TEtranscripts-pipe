import gzip
import itertools
import os

import pandas as pd
from snakemake.utils import validate

# How many FASTQ records to peek at when auto-detecting read length for
# STAR's --sjdbOverhang (see _fastq_read_length / _resolve_sjdb_overhang).
_FASTQ_READ_LENGTH_N_READS = 20

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

V = config["versions"]

# -----------------------------------------------------------------------------
# Per-rule compute resources (threads/mem_mb/runtime), loaded from
# workflow/default-config/resources.yaml (built-in defaults) and optionally
# overridden by config/resources.yaml if you create one. A rule (or a
# missing key within its entry) not
# present there falls back to a small conservative default instead of
# failing -- see the "HPC / SLURM" section of the README for how these feed
# into cluster execution.
# -----------------------------------------------------------------------------
RESOURCES = config.get("resources", {})
_RESOURCE_DEFAULTS = RESOURCES.get("_default", {"threads": 1, "mem_mb": 4000, "runtime": 60})


def get_resources(rule_name):
    """Return {threads, mem_mb, runtime} for a rule name."""
    return {**_RESOURCE_DEFAULTS, **RESOURCES.get(rule_name, {})}

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
# empty optional cells (fastq_2, strandedness, condition) need to be
# converted explicitly. astype(object) first is required on pandas>=2's
# str-backed columns, which otherwise silently keep NaN instead of None.
samples = samples.astype(object).where(pd.notnull(samples), None)
validate(samples, schema="../schemas/samples.schema.yaml")

SAMPLES = list(samples["sample"])
HAS_CONDITION = "condition" in samples.columns and samples["condition"].notna().all()

# Check sample-sheet paths exist at parse time so missing files fail fast.
_missing_paths = []
for col in ("fastq_1", "fastq_2"):
    if col in samples.columns:
        for sample, path in samples[col].items():
            if pd.notna(path) and str(path).strip():
                if not os.path.exists(str(path)):
                    _missing_paths.append((sample, col, path))
if _missing_paths:
    raise ValueError(
        "Sample sheet references files that do not exist:\n"
        + "\n".join(f"  {s}.{c}: {p}" for s, c, p in _missing_paths)
        + "\nRun Snakemake from the directory these paths are relative to "
        "(usually the repo root), or fix the paths in the sample sheet."
    )

# -----------------------------------------------------------------------------
# Reference files: transparently support gzipped fasta/gtf/te_gtf.
# STAR's --genomeFastaFiles/--sjdbGTFfile and TEcount/TEtranscripts'
# --GTF/--TE all expect plain-text files, so any ref.* path ending in .gz is
# decompressed once via the gunzip_reference rule (ref.smk); everything else
# in the workflow refers to the resolved (always-uncompressed) path below,
# never to config["ref"][...] directly.
# -----------------------------------------------------------------------------
REFERENCE_GZ_SOURCES = {}  # decompressed stem -> original .gz path, for gunzip_reference


def _resolve_ref_path(key):
    path = config["ref"][key]
    if str(path).endswith(".gz"):
        stem = os.path.basename(str(path))[:-3]  # strip trailing ".gz"
        decompressed = f"resources/decompressed/{stem}"
        if stem in REFERENCE_GZ_SOURCES and REFERENCE_GZ_SOURCES[stem] != path:
            raise ValueError(
                f"Two different ref.* .gz files decompress to the same "
                f"filename '{stem}' (resources/decompressed/{stem}) -- "
                "rename one of them so their basenames are unique."
            )
        REFERENCE_GZ_SOURCES[stem] = path
        return decompressed
    return path


FASTA = _resolve_ref_path("fasta")
GTF = _resolve_ref_path("gtf")
TE_GTF = _resolve_ref_path("te_gtf")

# -----------------------------------------------------------------------------
# STAR sjdbOverhang: auto-detected from the sample sheet's fastq files
# (max read length - 1, STAR's own recommendation) unless the user pins an
# explicit integer in config.yaml. Detection peeks at the first few reads of
# every fastq_1/fastq_2 file, gzip-aware, so it works whether reads are
# compressed or not.
# -----------------------------------------------------------------------------
def _open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _fastq_read_length(path, n_reads=_FASTQ_READ_LENGTH_N_READS):
    """Max read length in the first n_reads of a fastq(.gz) file, or 0 if
    the file doesn't exist / can't be read (e.g. not downloaded yet)."""
    if not path or not os.path.exists(path):
        return 0
    max_len = 0
    try:
        with _open_maybe_gz(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:  # sequence line of each fastq record
                    max_len = max(max_len, len(line.strip()))
                if i >= n_reads * 4:
                    break
    except OSError:
        return 0
    return max_len


def _resolve_sjdb_overhang():
    configured = config["ref"].get("sjdb_overhang", "auto")
    if str(configured).strip().lower() != "auto":
        return int(configured), "config.yaml (explicit)"
    paths = []
    for col in ("fastq_1", "fastq_2"):
        if col in samples.columns:
            paths.extend([p for p in samples[col].tolist() if p])
    max_len = max((_fastq_read_length(p) for p in paths), default=0)
    if max_len == 0:
        cwd = os.getcwd()
        checked = "\n".join(
            f"  - {p}  (exists: {os.path.exists(p)}, resolves to: "
            f"{os.path.abspath(p)})"
            for p in paths
        ) or "  (no fastq_1/fastq_2 values found in the sample sheet at all)"
        raise ValueError(
            "ref.sjdb_overhang is 'auto' but no readable fastq file was "
            f"found to auto-detect read length from. Current working "
            f"directory (paths are resolved relative to this): {cwd}\n"
            f"Paths checked from {config['samples']}:\n{checked}\n"
            "Either run Snakemake from the directory these paths are "
            "relative to (usually the repo root), fix the paths in the "
            "sample sheet, or set ref.sjdb_overhang to an explicit integer "
            "in config.yaml to skip auto-detection entirely."
        )
    return max_len - 1, f"auto-detected (max read length {max_len} - 1)"


SJDB_OVERHANG, _SJDB_OVERHANG_SOURCE = _resolve_sjdb_overhang()

# Record how key values were resolved -- a lightweight, always-current log of
# workflow parse-time decisions, separate from the per-rule logs in logs/.
os.makedirs("logs", exist_ok=True)
with open("logs/config_resolution.log", "w") as fh:
    fh.write(f"sjdb_overhang = {SJDB_OVERHANG} ({_SJDB_OVERHANG_SOURCE})\n")
    if REFERENCE_GZ_SOURCES:
        fh.write("gzipped reference files to be decompressed:\n")
        for stem, gz_path in REFERENCE_GZ_SOURCES.items():
            fh.write(f"  {gz_path} -> resources/decompressed/{stem}\n")
    else:
        fh.write("no gzipped reference files (fasta/gtf/te_gtf) detected\n")

# -----------------------------------------------------------------------------
# Per-sample strandedness resolution.
# The sample sheet's optional "strandedness" column (nf-core/rnaseq
# vocabulary: auto/forward/reverse/unstranded, so an nf-core samplesheet can
# be used as-is) takes priority per-sample; if empty/absent for a sample, it
# defaults to "auto". Resolved and validated once, eagerly, so a typo fails
# fast at parse time rather than mid-run.
# -----------------------------------------------------------------------------
STRANDEDNESS_ALIASES = {
    "auto": "auto",
    "forward": "forward",
    "reverse": "reverse",
    "unstranded": "no",  # nf-core vocabulary -> TEtranscripts vocabulary
    "no": "no",  # TEtranscripts' own vocabulary, accepted directly too
}


def _normalize_stranded(raw, context=""):
    key = str(raw).strip().lower()
    if key not in STRANDEDNESS_ALIASES:
        raise ValueError(
            f"Unrecognized strandedness value '{raw}'{context}. Expected one "
            "of: auto, forward, reverse, unstranded (nf-core vocabulary), or "
            "no (TEtranscripts vocabulary)."
        )
    return STRANDEDNESS_ALIASES[key]


# Default when a sample sheet has no "strandedness" column, or leaves it
# blank for a given sample: auto-detect via RSeQC.
DEFAULT_STRANDEDNESS_MODE = "auto"


def _sample_effective_mode(sample):
    if "strandedness" in samples.columns:
        val = samples.loc[sample, "strandedness"]
        if pd.notna(val) and str(val).strip() != "":
            return _normalize_stranded(val, f" for sample '{sample}' in {config['samples']}")
    return DEFAULT_STRANDEDNESS_MODE


SAMPLE_STRANDED_MODE = {s: _sample_effective_mode(s) for s in SAMPLES}
AUTO_SAMPLES = [s for s in SAMPLES if SAMPLE_STRANDED_MODE[s] == "auto"]


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
    """--readFilesCommand zcat if both fastq files are gzipped."""
    d = star_input(wildcards)
    fq1_gz = str(d["fq1"]).endswith(".gz")
    fq2_gz = "fq2" not in d or str(d["fq2"]).endswith(".gz")
    if fq1_gz and fq2_gz:
        return "--readFilesCommand zcat"
    if fq1_gz != fq2_gz:
        raise ValueError(
            f"Mixed compression: fastq_1 is gzipped but fastq_2 is not "
            f"(or vice versa) for sample '{wildcards.sample}'. "
            "Both files must be either both gzipped or both uncompressed."
        )
    return ""


def strandedness_input(wildcards):
    """Dependency on this sample's auto-detected strandedness call.

    Returns [] (no dependency) when this sample's effective mode is fixed
    (from the sample sheet's "strandedness" column or the config default),
    since no RSeQC/auto-detection step is needed in that case.
    """
    if SAMPLE_STRANDED_MODE[wildcards.sample] == "auto":
        return f"results/rseqc/{wildcards.sample}/strandedness.txt"
    return []


def get_strandedness_param(wildcards, input):
    """Resolve the --stranded value (no/forward/reverse) for a sample.

    If this sample's effective mode is "auto", read the value that
    workflow/scripts/determine_strandedness.py determined from the RSeQC
    infer_experiment.py output (the file is guaranteed to already exist
    because it is declared as a rule input alongside this params function).
    Otherwise, the fixed value (sample sheet override or config default) is
    used directly.
    """
    mode = SAMPLE_STRANDED_MODE[wildcards.sample]
    if mode == "auto":
        with open(input.strandedness) as fh:
            value = fh.read().strip()
        if value not in ("no", "forward", "reverse"):
            raise ValueError(
                f"Invalid strandedness value '{value}' in "
                f"{input.strandedness} for sample '{wildcards.sample}'. "
                "Expected one of: no, forward, reverse."
            )
        return value
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


def _contrast_samples(wildcards):
    c = CONTRASTS[wildcards.contrast]
    return list(c["treatment"]) + list(c["control"])


def contrast_strandedness_input(wildcards):
    """Dependency on the auto-detection call of every sample in a contrast
    whose effective strandedness mode is "auto" (others are fixed and need
    no RSeQC run)."""
    auto_samples = [s for s in _contrast_samples(wildcards) if SAMPLE_STRANDED_MODE[s] == "auto"]
    return expand("results/rseqc/{sample}/strandedness.txt", sample=auto_samples)


def get_contrast_strandedness_param(wildcards, input):
    """Resolve a single --stranded value shared by every sample in a contrast.

    TEtranscripts runs one DESeq2 analysis across all treatment/control BAMs
    at once, so it needs one strandedness value. Each sample's own effective
    mode (sample-sheet override, auto-detected, or config default) is
    resolved and, if they disagree, this fails loudly rather than silently
    picking one -- that mismatch usually means samples were prepared with
    different library kits and shouldn't be pooled blindly.
    """
    all_samples = _contrast_samples(wildcards)
    auto_samples = [s for s in all_samples if SAMPLE_STRANDED_MODE[s] == "auto"]
    file_map = dict(zip(auto_samples, input.strandedness))
    values = set()
    for s in all_samples:
        mode = SAMPLE_STRANDED_MODE[s]
        if mode == "auto":
            with open(file_map[s]) as fh:
                value = fh.read().strip()
            if value not in ("no", "forward", "reverse"):
                raise ValueError(
                    f"Invalid strandedness value '{value}' in "
                    f"{file_map[s]} for sample '{s}'. "
                    "Expected one of: no, forward, reverse."
                )
            values.add(value)
        else:
            values.add(mode)
    if len(values) > 1:
        raise ValueError(
            f"Samples in contrast '{wildcards.contrast}' have inconsistent "
            f"strandedness ({values}). Check the 'strandedness' column in "
            f"{config['samples']}, or verify auto-detected values agree."
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
# Container image used for every rule (replaces per-rule conda environments).
# When the workflow is run with --sdm apptainer / --sdm singularity / --sdm
# docker, Snakemake executes every shell command inside this container.
#
# Default / local development: build from the Dockerfile.
# Production/HPC: pull from a pre-built image on Docker Hub (faster).
# Switch between them here or override CONTAINER_IMAGE in a rule.
# -----------------------------------------------------------------------------
CONTAINER_IMAGE = "docker://altintasali/rnaseq-star-tetranscripts:latest"

container:
    CONTAINER_IMAGE
