# Contrasts, the all_*() target lists, benchmark bookkeeping.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

# -----------------------------------------------------------------------------
# Pairwise contrasts, derived automatically from the "condition" column.
# If the sample sheet has no "condition" column, CONTRASTS is empty and the
# TEtranscripts differential-analysis step is skipped entirely.
# -----------------------------------------------------------------------------
def _build_contrasts():
    if not HAS_CONDITION:
        return {}
    conditions = sorted(samples["condition"].dropna().unique())
    contrasts = {}
    for control, treatment in itertools.combinations(conditions, 2):
        name = f"{treatment}_vs_{control}"
        contrasts[name] = {
            "treatment": list(samples.loc[samples["condition"] == treatment, "sample"]),
            "control": list(samples.loc[samples["condition"] == control, "sample"]),
        }
    return contrasts


CONTRASTS = _build_contrasts()

# Warn at parse time when any contrast has fewer than 2 replicates in a
# group -- DESeq2 (called internally by TEtranscripts) needs >= 2
# replicates per condition to estimate dispersion.  This is a warning,
# not an error, because the user may be doing exploratory analysis.
import warnings as _warnings

for _cname, _cval in CONTRASTS.items():
    for _group in ("treatment", "control"):
        if len(_cval[_group]) < 2:
            _warnings.warn(
                f"contrast '{_cname}': {_group} group has only "
                f"{len(_cval[_group])} sample(s); DESeq2 needs >= 2 "
                f"replicates per condition for dispersion estimation.",
                stacklevel=2,
            )


def _contrast_samples(wildcards):
    c = CONTRASTS[wildcards.contrast]
    return list(c["treatment"]) + list(c["control"])


def contrast_strandedness_input(wildcards):
    """Dependency on the auto-detection call of every sample in a contrast
    whose effective strandedness mode is "auto" (others are fixed and need
    no RSeQC run)."""
    auto_samples = [s for s in _contrast_samples(wildcards) if SAMPLE_STRANDED_MODE[s] == "auto"]
    return expand("results/rseqc/{sample}_strandedness.txt", sample=auto_samples)


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
                values.add(fh.read().strip())
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
    return expand("results/tecount/{sample}.cntTable.gz", sample=SAMPLES)


def all_telocal_outputs():
    """TElocal per-sample count tables, counts matrix + sample-QC artifacts,
    and summary barplots for the `all` target (Snakefile). The QC view only
    runs when telocal.qc.enabled is true; the summary barplots always
    render."""
    files = expand("results/telocal/{sample}.cntTable.gz", sample=SAMPLES)
    if TELOCAL_QC_ENABLED:
        transform = TELOCAL_QC["pca_transform"]
        files += [
            "results/telocal/counts_matrix.tsv.gz",
            f"results/telocal/qc/{transform}_counts.tsv.gz",
            f"results/telocal/qc/pca_{transform}_mqc.json",
            f"results/telocal/qc/heatmap_{transform}_mqc.json",
        ]
    files += [
        "results/telocal/qc/telocal_assignment_mqc.json",
        "results/telocal/qc/telocal_te_class_mqc.json",
        "results/telocal/telocal_locations.bed",
    ]
    return files


def all_trim_outputs():
    """One trimmed fastq path per sample (only when trimming is enabled),
    used by the `trimming_only` convenience target (Snakefile)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1.fq.gz"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed.fq.gz"
        )
        for s in SAMPLES
    ]


def all_fastqc_reports():
    """One FastQC .zip report path per trimmed sample (only when trimming is
    enabled), used by the multiqc rule so it scans results/trimming/ for the
    FastQC + TrimGalore! reports without depending on the trimmed fastq files
    themselves (which may be temp()-deleted after alignment)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1_fastqc.zip"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed_fastqc.zip"
        )
        for s in SAMPLES
    ]


def all_raw_fastqc_reports():
    """FastQC .zip report paths for the raw (merged) input fastqs, one per
    sample/read. Run unconditionally so the MultiQC report always covers the
    untrimmed input, regardless of the optional `trimming` step."""
    return [
        f"results/fastqc/raw/{s}_R{read}_fastqc.zip"
        for s in SAMPLES
        for read in (1, 2)
        if read == 1 or _is_paired(s)
    ]


def all_diffexp_outputs():
    return expand(
        "results/tetranscripts/{contrast}_sigdiff_gene_TE.txt.gz",
        contrast=list(CONTRASTS.keys()),
    )


def all_benchmark_files():
    """Every benchmark file this configuration will produce, so the
    benchmark_summary rule (qc.smk) aggregates exactly the rules that ran
    into the MultiQC resource-usage section. Kept in lockstep with the
    `benchmark:` declarations in the rules -- only the conditional ones
    (merging, trimming, gunzip, RSeQC bed12 conversion, strandedness
    auto-detection, contrasts) need to be gated here. The multiqc and
    benchmark_summary rules' own benchmarks are deliberately excluded to
    avoid a cyclic dependency (their resource use is negligible)."""
    files = [
        "results/pipeline_info/benchmarks/software_versions/software_versions.txt",
        "results/pipeline_info/benchmarks/gene_name_lookup/gene_name_lookup.txt",
    ]
    # The star_index rule only exists when we build the index ourselves; with
    # star.build_index: false the external index is a plain input and no such
    # job (or benchmark) is ever produced.
    if STAR_BUILD_INDEX:
        files.append(
            "results/pipeline_info/benchmarks/star_index/star_index.txt"
        )
    # BED12 gene-model conversion only runs when RSeQC strandedness
    # auto-detection actually needs it.
    if STRAND_CHECK_SAMPLES:
        files.append(
            "results/pipeline_info/benchmarks/strandedness_check/"
            "strandedness_check.txt"
        )
    if AUTO_SAMPLES:
        files += [
            "results/pipeline_info/benchmarks/gtf_to_genepred/gtf_to_genepred.txt",
            "results/pipeline_info/benchmarks/genepred_to_bed12/genepred_to_bed12.txt",
        ]
    for stem in REFERENCE_GZ_SOURCES:
        files.append(f"results/pipeline_info/benchmarks/gunzip_reference/{stem}.txt")
    for s in SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/star_align/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_sort/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_index/{s}.txt",
            f"results/pipeline_info/benchmarks/tecount/{s}.txt",
            f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R1.txt",
        ]
        if _is_paired(s):
            files.append(f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R2.txt")
        # Lane-merging (cat_fastq) only for samples with multiple lanes/runs.
        for read in (1, 2):
            if (read == 1 or _is_paired(s)) and len(sample_fastqs(s, read)) > 1:
                files.append(
                    f"results/pipeline_info/benchmarks/cat_fastq/{s}_R{read}.txt"
                )
        if TRIM_ENABLED:
            files.append(
                f"results/pipeline_info/benchmarks/"
                f"{'trim_galore_pe' if _is_paired(s) else 'trim_galore_se'}/{s}.txt"
            )
    for s in STRAND_CHECK_SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/rseqc_infer_experiment/{s}.txt",
            f"results/pipeline_info/benchmarks/determine_strandedness/{s}.txt",
        ]
    for contrast in CONTRASTS:
        files.append(
            f"results/pipeline_info/benchmarks/tetranscripts_diffexp/{contrast}.txt"
        )
    # Chimera-screen rules only run when the chimera stage is enabled.
    if CHIMERA_JUNCTION_ENABLED:
        files += [
            "results/pipeline_info/benchmarks/annotation_to_bed/annotation_to_bed.txt",
            "results/pipeline_info/benchmarks/chimera_counts/chimera_counts.txt",
            "results/pipeline_info/benchmarks/junction_highlights/"
            "junction_highlights.txt",
            "results/pipeline_info/benchmarks/chimera_evidence/"
            "chimera_evidence.txt",
            "results/pipeline_info/benchmarks/chimera_evidence_guide/"
            "chimera_evidence_guide.txt",
            "results/pipeline_info/benchmarks/chimera_candidates_table/"
            "chimera_candidates_table.txt",
            "results/pipeline_info/benchmarks/chimera_evidence_heatmap/"
            "chimera_evidence_heatmap.txt",
            "results/pipeline_info/benchmarks/sample_evidence_status/"
            "sample_evidence_status.txt",
        ]
        for s in SAMPLES:
            files += [
                f"results/pipeline_info/benchmarks/parse_chimeric_junctions/{s}.txt",
                f"results/pipeline_info/benchmarks/junction_qc/{s}.txt",
            ]
            if config["chimera"]["junction"]["outputs"]["write_igv_bed"]:
                files.append(
                    f"results/pipeline_info/benchmarks/chimera_igv_bed/{s}.txt"
                )
        if config["chimera"]["junction"]["outputs"]["write_counts_matrix"]:
            transform = config["chimera"]["junction"]["qc"]["pca_transform"]
            files += [
                f"results/pipeline_info/benchmarks/"
                f"sample_qc_transform/{transform}.txt",
                f"results/pipeline_info/benchmarks/sample_qc/{transform}.txt",
            ]
    # TEcounts sample-QC rules only run when tetranscripts.qc.enabled.
    if TECOUNT_QC_ENABLED:
        transform = TECOUNT_QC["pca_transform"]
        files += [
            "results/pipeline_info/benchmarks/tecount_counts/tecount_counts.txt",
            f"results/pipeline_info/benchmarks/tecount_qc_transform/{transform}.txt",
            f"results/pipeline_info/benchmarks/tecount_qc/{transform}.txt",
        ]
    # The summary-barplot rule runs on every run (raw cntTables only).
    files.append(
        "results/pipeline_info/benchmarks/tecount_summary/tecount_summary.txt"
    )
    # TElocal rules only run when the telocal stage is enabled; its QC view
    # additionally requires telocal.qc.enabled.
    if TELOCAL_ENABLED:
        files += [
            "results/pipeline_info/benchmarks/telocal_locind/locind.txt",
            "results/pipeline_info/benchmarks/telocal_summary/telocal_summary.txt",
            "results/pipeline_info/benchmarks/telocal_locations/locations.txt",
        ]
        for s in SAMPLES:
            files.append(f"results/pipeline_info/benchmarks/telocal/{s}.txt")
        if TELOCAL_QC_ENABLED:
            transform = TELOCAL_QC["pca_transform"]
            files += [
                "results/pipeline_info/benchmarks/telocal_counts/telocal_counts.txt",
                f"results/pipeline_info/benchmarks/telocal_qc_transform/{transform}.txt",
                f"results/pipeline_info/benchmarks/telocal_qc/{transform}.txt",
            ]
    return sorted(set(files))


def allocated_resources_by_rule():
    """{rule: {"threads", "mem_mb"}} for every rule that has benchmark files
    (the benchmark_summary rule's input), read from resources.yaml -- the
    per-job allocation against which the benchmark_summary script computes
    CPU/RAM efficiency. Iterated in Snakemake's own rule order (workflow.rules,
    an OrderedDict) so the resource-usage table lists rules in the order they
    appear in the workflow, not alphabetically."""
    benchmark_rules = {
        os.path.basename(os.path.dirname(path)) for path in all_benchmark_files()
    }
    out = {}
    for r in workflow.rules:
        if r.name in benchmark_rules:
            out[r.name] = get_resources(r.name)
    return out
