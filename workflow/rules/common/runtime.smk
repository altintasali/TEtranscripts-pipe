# Star 2-pass, strandedness resolution, fastq preparation.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

# -----------------------------------------------------------------------------
# STAR 2-pass mapping (STAR manual section 9): "none" (single-pass),
# "per_sample" (--twopassMode Basic, no other rules change), or "cohort"
# (novel junctions pooled across all samples before any sample's final
# alignment -- see rules/star_two_pass.smk, only included when this is
# "cohort"). Schema enum already rejects any other value.
#
# Default is "cohort" FOR NOW while this feature is being evaluated -- set
# star.two_pass: none explicitly to opt out and get the original,
# single-pass-only DAG.
# -----------------------------------------------------------------------------
STAR_TWO_PASS = config["star"].get("two_pass", "cohort")
if STAR_TWO_PASS in ("per_sample", "cohort"):
    logger.warning(
        f"star.two_pass is {STAR_TWO_PASS!r}. STAR 2-pass mapping roughly "
        "doubles STAR's alignment runtime per sample -- if star_align jobs "
        "start timing out, increase its runtime in input/resources.yaml."
        + (
            " Cohort mode also adds star_align_pass1/star_merge_junctions "
            "jobs and a DAG sync point: no sample's final alignment starts "
            "until every sample's pass-1 alignment has finished."
            if STAR_TWO_PASS == "cohort"
            else ""
        )
    )

if STAR_BUILD_INDEX:
    # We build it, so the workflow owns the directory: it is a rule OUTPUT.
    STAR_INDEX_DIR = os.path.abspath(_STAR_INDEX_RAW or "results/star_index")
else:
    # External, pre-built index: used directly, as a rule INPUT.
    #
    # It is deliberately NOT a rule output. Snakemake deletes a rule's
    # directory() output before running the job, so while star_index
    # declared this path as its output, pointing star.index at a real index
    # meant Snakemake owned -- and wiped -- it. The previous workaround was
    # to copy the whole index into results/ so the rule's output always sat
    # inside the workflow; that cost a full copy of a multi-GB index on
    # every run, and made safety positional (get the variable wrong and it
    # is rm -rf again) rather than structural.
    #
    # Now star_index is simply not defined when build_index is false (see
    # ref.smk), so no rule produces this path and Snakemake treats it as a
    # required pre-existing input. Inputs are never deleted, so the original
    # cannot be touched by any code path -- and there is no copy.
    STAR_INDEX_DIR = os.path.abspath(_STAR_INDEX_RAW) if _STAR_INDEX_RAW else ""

    if not STAR_INDEX_DIR:
        raise WorkflowError(
            "star.build_index is false but no star.index path was provided. "
            "Set star.index in config.yaml to the directory containing the "
            "pre-built STAR index."
        )
    if not os.path.isdir(STAR_INDEX_DIR):
        raise WorkflowError(
            f"star.build_index is false but star.index path does not exist: "
            f"{STAR_INDEX_DIR}. Either build the index beforehand or set "
            f"star.build_index: true in config.yaml."
        )
    # STAR needs these three to load a genome; catching it here beats a
    # cryptic STAR failure once jobs are already queued.
    _missing = [
        f for f in ("SA", "SAindex", "Genome")
        if not os.path.isfile(os.path.join(STAR_INDEX_DIR, f))
    ]
    if _missing:
        raise WorkflowError(
            f"star.index ({STAR_INDEX_DIR}) does not look like a STAR index: "
            f"missing {', '.join(_missing)}. Point star.index at the "
            f"directory produced by STAR --runMode genomeGenerate, or set "
            f"star.build_index: true to build one."
        )

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

# Samples RSeQC actually runs on. AUTO_SAMPLES are the ones that NEED it --
# their --stranded value comes from the inference. With check_provided (the
# default) it runs on the rest too, purely so the report can compare a
# declared strandedness against the data and flag a disagreement: an
# explicitly set sample sheet value is never overridden by the inference
# (see get_strandedness_param, which dispatches on SAMPLE_STRANDED_MODE and
# does not look at the RSeQC call for a non-auto sample).
#
# A wrong strandedness silently halves or inverts every count downstream
# without failing anything, which is exactly the class of error worth
# spending a cheap extra job to catch.
STRAND_CHECK_PROVIDED = bool(
    config.get("strandedness", {}).get("check_provided", True)
)
STRAND_CHECK_SAMPLES = SAMPLES if STRAND_CHECK_PROVIDED else AUTO_SAMPLES


def _is_paired(sample):
    """A sample is paired-end if it has at least one non-empty fastq_2 value.

    This is resolved per-sample (not globally), so a single sample sheet can
    freely mix paired-end and single-end samples. Every lane/run of one
    sample must agree (enforced while loading the sample sheet above).
    """
    return len(sample_fastqs(sample, 2)) > 0


# -----------------------------------------------------------------------------
# Fastq preparation pipeline: raw lanes -> (merge) -> (trim) -> STAR.
# Samples sequenced across multiple lanes/runs are concatenated by
# rules/fastq.smk (cat_fastq) into results/fastq/. If trimming is enabled,
# rules/trimming.smk (trim_galore, nf-core/rnaseq defaults) then runs on the
# merged files and STAR consumes the trimmed output; otherwise STAR reads the
# merged (or raw, for single-lane) files directly.
# -----------------------------------------------------------------------------
TRIM_ENABLED = bool(config.get("trimming", {}).get("enabled", True))

# Optional chimera-junction screen (rules/chimera_reads.smk +
# chimera_reads_qc.smk). When enabled (the default), the STAR alignment emits
# chimeric junctions and the chimera rules annotate them; set
# chimera.reads.enabled: false to opt out -- no chimera STAR flags are
# passed, the chimera rules are not included, and the workflow behaves like
# the plain quantification pipeline. See CHIMERA_ASSEMBLY_ENABLED below for
# the complementary StringTie-assembly-based screen.
CHIMERA_READS_ENABLED = bool(
    config.get("chimera", {}).get("reads", {}).get("enabled", True)
)

# Optional chimera-assembly screen (rules/chimera_assembly.smk):
# StringTie-assembly-based detection, complementing CHIMERA_READS_ENABLED
# above. Off by default -- newer and less validated.
CHIMERA_ASSEMBLY_ENABLED = bool(
    config.get("chimera", {}).get("assembly", {}).get("enabled", False)
)

# Sample-QC thresholds for the chimera views (PCA / sample clustering). Lives
# here, not in chimera_reads.smk, because BOTH screens' QC views use it and
# that file is included only when the junction screen is on -- referencing it
# from chimera_assembly.smk would NameError on an assembly-only run, which is
# exactly the configuration guard 27 pins.
#
# The assembly view borrows these rather than having a parallel config block:
# the two views answer the same question and there is no evidence they want
# different cut-offs. Split them if that stops being true.
CHIMERA_QC = config["chimera"]["reads"]["qc"]

# How many candidate rows the report's table renders. MultiQC embeds table
# data in the HTML and a real cohort produces tens of thousands of gene-TE
# pairs, so the section shows a head and points at candidates.tsv.gz for the
# rest. Not a config key: a rendering limit, not an analysis choice.
CHIMERA_TABLE_TOP_N = 50

# TEcounts sample-QC (PCA + sample clustering, rules/tecount_qc.smk), built
# from the per-sample TEcount tables. Defaults come from the built-in
# workflow/default-config/tetranscripts.yaml (the `qc:` section); user config
# overrides them. When disabled, no counts matrix or QC plots are produced.
TECOUNT_QC = config["tetranscripts"]["qc"]
TECOUNT_QC_ENABLED = bool(TECOUNT_QC["enabled"])

# TElocal locus-level quantification (rules/telocal.smk). Requires a
# pre-built .locInd file; disabled by default in older configs but enabled
# in the built-in telocal.yaml defaults.  When locind is empty, the
# telocal_locind rule auto-builds from the TE GTF.
TELOCAL_ENABLED = bool(config.get("telocal", {}).get("enabled", False))
_telocal_locind_cfg = config.get("telocal", {}).get("locind", "")
TELOCAL_QC = config["telocal"].get(
    "qc",
    {
        "enabled": True,
        "min_samples_present": 2,
        "min_total_counts": 5,
        "min_events": 10,
        "pca_transform": "log2",
        "feature_class": "TE",
    },
)
TELOCAL_QC_ENABLED = bool(TELOCAL_QC["enabled"])

if TELOCAL_ENABLED and TELOCAL_QC_ENABLED and TELOCAL_QC["pca_transform"] in ("vst", "rlog"):
    logger.warning(
        f"telocal.qc.pca_transform is {TELOCAL_QC['pca_transform']!r}. DESeq2's "
        "vst/rlog on TElocal's locus-level matrix (one row per TE instance, "
        "often millions) can be prohibitively slow -- the default is 'log2' "
        "for this reason (see README). telocal_qc_transform's runtime "
        "budget (workflow/default-config/resources.yaml, or your "
        "input/resources.yaml override) may need increasing well beyond its "
        "default if the job times out."
    )


def _telocal_locind_path():
    """Return the .locInd path for the telocal rule.

    When the user provides a path, use it directly.  When locind is empty,
    auto-build from the TE GTF (telocal_locind rule) into results/telocal/.
    The name must keep the .locInd suffix -- TElocal rejects any --TE file
    whose path does not end in .locInd.
    """
    if _telocal_locind_cfg:
        return _telocal_locind_cfg
    return "results/telocal/telocal.locInd"

# TrimGalore! always appends _trimmed (single-end) or _val_1/_val_2 (paired)
# to the *input* basename, and its --basename normalization only strips a
# single "_trimmed"/"_val_1" suffix. Feeding it an already-trimmed fastq
# therefore produces a no-op rename like foo_trimmed_trimmed_trimmed.fq.gz
# while the declared results/trimming/<sample>_trimmed.fq.gz output is never
# created (or a stale file is silently reused). Fail fast at parse time
# instead of surfacing a confusing MissingOutputException mid-run.
if TRIM_ENABLED:
    _TRIMMED_FASTQ_RE = re.compile(
        r"(_trimmed|_val_[12]|\.trimmed)(\.fq|\.fastq)(\.gz)?$", re.IGNORECASE
    )
    _offenders = [
        f"{sample}: {path}"
        for sample in SAMPLES
        for read in (1, 2)
        for path in sample_fastqs(sample, read)
        if _TRIMMED_FASTQ_RE.search(os.path.basename(str(path)))
    ]
    if _offenders:
        raise WorkflowError(
            "Trimming is enabled (trimming.enabled), but the sample sheet "
            "points at fastqs that already look trimmed:\n  "
            + "\n  ".join(_offenders)
            + "\nFeeding already-trimmed reads to TrimGalore! produces a "
            "no-op rename and breaks the pipeline. Point the sample sheet at "
            "untrimmed reads, or set `trimming.enabled: false` in config.yaml "
            "to use these files as-is."
        )

# Whether the merged (lane-concatenated) and trimmed fastq files are kept
# after downstream rules (STAR/alignment) are done with them. Set either to
# false in the `outputs` config section and those files are marked temp() --
# Snakemake deletes them after the last consumer finishes, saving disk on big
# runs at the cost of re-running the merge/trim steps if they're needed again.
KEEP_MERGED_FASTQ = bool(config.get("outputs", {}).get("keep_merged_fastq", True))
KEEP_TRIMMED_FASTQ = bool(config.get("outputs", {}).get("keep_trimmed_fastq", True))
KEEP_STAR_INDEX = bool(config.get("outputs", {}).get("keep_star_index", True))
# Whether the auto-built TElocal .locInd index
# (results/telocal/telocal.locInd) is
# kept after all TElocal runs finish. Never applies to a user-provided
# telocal.locind path -- that file is only ever read, never deleted.
KEEP_TELOCAL_INDEX = bool(
    config.get("outputs", {}).get("keep_telocal_index", True)
)


def _maybe_temp(path, keep):
    """A rule output path, wrapped in temp() unless the user opted to keep
    it via the `outputs` config section. Snakemake allows temp() on wildcard
    patterns, so this is applied at parse time (the keep flag is a constant)."""
    if keep:
        return path
    return temp(path)


# The cat_fastq rule's output pattern, temp()-wrapped when the merged fastqs
# are not kept (see the `outputs` config section).
MERGED_FASTQ_OUTPUT = _maybe_temp(
    "results/fastq/{sample}_R{read}.fastq.gz", KEEP_MERGED_FASTQ
)


# Directory for STAR's per-run scratch files (the _STARtmp directory), outside
# the results tree. See the star_align rule in align.smk.
STAR_TMPDIR = config["star"].get("tmpdir", "/tmp")


def merged_fastq_path(sample, read):
    """Path of the concatenated fastq for a sample's read (results/fastq/).

    Samples with a single lane/run have no merge step: their raw path is
    returned directly so no unnecessary copy is made. Multi-lane samples are
    always merged to a .gz file (plain-text lanes are gzipped on the fly by
    the cat_fastq rule), so STAR's --readFilesCommand zcat applies uniformly.
    """
    files = sample_fastqs(sample, read)
    if len(files) <= 1:
        return files[0] if files else None
    return f"results/fastq/{sample}_R{read}.fastq.gz"


def star_fastq(sample, read):
    """The fastq path STAR aligns for a sample's read (after merge/trim)."""
    base = merged_fastq_path(sample, read)
    if base is None:
        return None
    if TRIM_ENABLED:
        if _is_paired(sample):
            return f"results/trimming/{sample}_val_{read}.fq.gz"
        return f"results/trimming/{sample}_trimmed.fq.gz"
    return base


def star_input(wildcards):
    """Return {fq1: ..., [fq2: ...]} paths for a sample."""
    d = {"fq1": star_fastq(wildcards.sample, 1)}
    if _is_paired(wildcards.sample):
        d["fq2"] = star_fastq(wildcards.sample, 2)
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
    """Dependency on this sample's auto-detected strandedness call.

    Returns [] (no dependency) when this sample's effective mode is fixed
    (from the sample sheet's "strandedness" column or the config default),
    since no RSeQC/auto-detection step is needed in that case.
    """
    if SAMPLE_STRANDED_MODE[wildcards.sample] == "auto":
        return f"results/rseqc/{wildcards.sample}_strandedness.txt"
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
            return fh.read().strip()
    return mode
