# Sample sheet load/validate, samples, per-sample fastqs.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

# -----------------------------------------------------------------------------
# Load & validate sample sheet
# columns: sample, fastq_1, fastq_2 (optional), strandedness (optional),
#          condition (optional).
#
# Repeated rows with the same `sample` name are allowed (nf-core/rnaseq
# style): they mean that sample's reads are split across lanes/runs and will
# be concatenated into one file per read before alignment (rules/fastq.smk).
# -----------------------------------------------------------------------------
if not os.path.exists(config["samples"]):
    raise WorkflowError(
        f"samples: '{config['samples']}' (from config.yaml, or --configfile/"
        "--config on the command line) does not exist.\n"
        f"Current working directory (relative paths are resolved against "
        f"this): {os.getcwd()}\n"
        f"That path resolves to: {os.path.abspath(config['samples'])}\n"
        "Snakemake resolves relative paths against the directory you *run "
        "the command from*, not the directory workflow/Snakefile lives in "
        "-- cd into the repo root (the directory that directly contains "
        "workflow/, config/, and .tests/) and try again."
    )


def _nonempty(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text or None


raw_samples = pd.read_csv(config["samples"], dtype=str, comment="#")
# jsonschema's "null" type only matches Python None, not pandas/NumPy NaN, so
# empty optional cells (fastq_2, strandedness, condition) need to be
# converted explicitly. astype(object) first is required on pandas>=2's
# str-backed columns, which otherwise silently keep NaN instead of None.
# Note: validate() below re-introduces NaN for empty cells (observed with
# snakemake 9.x), so _nonempty() also treats NaN/NA as empty defensively.
raw_samples = raw_samples.astype(object).where(pd.notnull(raw_samples), None)
validate(raw_samples, schema="../../schemas/samples.schema.yaml")

# Group repeated rows for the same sample into one per-sample entry holding a
# list of fastq paths per read. Consistency is enforced eagerly at parse time
# so a typo in a lane-split sample fails fast instead of mid-run.
_rows = {}
for _, row in raw_samples.iterrows():
    sample = _nonempty(row["sample"])
    if sample is None:
        raise WorkflowError(f"{config['samples']}: a row has no sample name.")
    entry = _rows.setdefault(
        sample,
        {"fastq_1": [], "fastq_2": [], "strandedness": None, "condition": None},
    )
    fq1 = _nonempty(row["fastq_1"])
    if fq1 is None:
        raise WorkflowError(f"sample '{sample}' has a row with no fastq_1.")
    entry["fastq_1"].append(fq1)
    fq2 = _nonempty(row["fastq_2"])
    if fq2 is not None:
        entry["fastq_2"].append(fq2)
    for col in ("strandedness", "condition"):
        value = _nonempty(row.get(col))
        if value is None:
            continue
        if entry[col] is not None and entry[col] != value:
            raise WorkflowError(
                f"sample '{sample}' has inconsistent '{col}' values across "
                f"rows ({entry[col]!r} vs {value!r}) in {config['samples']}."
            )
        entry[col] = value

for sample, entry in _rows.items():
    n_r1, n_r2 = len(entry["fastq_1"]), len(entry["fastq_2"])
    if n_r2 and n_r2 != n_r1:
        raise WorkflowError(
            f"sample '{sample}' has {n_r1} fastq_1 files but {n_r2} fastq_2 "
            f"files in {config['samples']} -- every lane/run must provide "
            "both reads (or the sample must be single-end throughout)."
        )

samples = (
    pd.DataFrame([{"sample": s, **e} for s, e in _rows.items()])
    .set_index("sample", drop=False)
    .sort_index()
)

SAMPLES = list(samples["sample"])
HAS_CONDITION = "condition" in samples.columns and samples["condition"].notna().all()


def sample_fastqs(sample, read):
    """Raw fastq paths for a sample's read 1/2, across all lanes/runs."""
    files = samples.loc[sample, f"fastq_{read}"]
    if isinstance(files, list):
        return [f for f in files if f]
    return [files] if files else []

# Warn (not fail) about sample-sheet fastq paths that don't exist yet: a typo
# would otherwise only surface when that rule's job runs -- potentially hours
# into a SLURM queue. Warning keeps pre-download dry-runs and laptops running
# against cluster-only paths working (the sjdb_overhang "auto" check above
# similarly tolerates missing files if an explicit value is given).
_missing_fastqs = [
    f"{sample}: {path}"
    for sample in SAMPLES
    for read in (1, 2)
    for path in sample_fastqs(sample, read)
    if path and not os.path.exists(path)
]
if _missing_fastqs:
    logger.warning(
        f"{len(_missing_fastqs)} fastq file(s) referenced by the sample sheet "
        "do not exist yet:\n  "
        + "\n  ".join(_missing_fastqs)
        + "\nThey must be present before the workflow can run; check the "
        "paths (and that you're in the directory they're relative to)."
    )
