def _fastqc_base_name(fname):
    """The name FastQC derives from an input fastq for its output .zip:
    the basename with recognized file extensions stripped (a leading .gz
    and then .fastq/.fq/.fasta/.fa/.txt). Used to find the zip FastQC
    produced so it can be renamed to the canonical {sample}_R{read} name."""
    base = os.path.basename(fname)
    for ext in (".gz", ".bz2", ".fastq", ".fq", ".fasta", ".fa", ".txt"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


rule fastqc_raw:
    # Always-on FastQC of the raw (merged) input fastqs, independent of the
    # optional `trimming` step, so the MultiQC report covers the untrimmed
    # input even when trimming is disabled.
    #
    # FastQC bakes the INPUT file's own basename into its report content --
    # the "Filename" line inside fastqc_data.txt, packed into the zip -- and
    # MultiQC's fastqc module unconditionally reads THAT to name the sample
    # in General Stats, not the zip's own filename. merged_fastq_path()
    # deliberately returns a single-lane sample's raw sample-sheet path
    # unchanged (no unnecessary copy), which can have an arbitrary basename
    # (e.g. from the sequencer), so without this, FastQC's embedded name --
    # and therefore the MultiQC sample name -- would not match the
    # {sample}_R{read} every other tool uses, breaking General Stats
    # alignment. Renaming just the output .zip (still done below, now a
    # no-op in the normal case) never touched this. Symlinking the input to
    # a canonically-named path before running FastQC fixes it at the
    # source -- verified empirically that FastQC's embedded Filename is the
    # symlink's own name, not its target's.
    input:
        fq=lambda wc: merged_fastq_path(wc.sample, int(wc.read)),
    output:
        zip="results/fastqc/raw/{sample}_R{read}_fastqc.zip",
    params:
        outdir="results/fastqc/raw",
        # The symlink FastQC actually reads, named the way every other tool
        # names this sample. Its extension is carried over from the real input
        # so FastQC still recognises the format.
        canonical=lambda wc, input: os.path.join(
            "results/fastqc/raw",
            f"{wc.sample}_R{wc.read}"
            + os.path.basename(input.fq)[len(_fastqc_base_name(input.fq)):],
        ),
        # Resolved here rather than in the shell: `realpath` is not portable
        # and the input may be a relative sample-sheet path.
        fq_abs=lambda wc, input: os.path.abspath(input.fq),
    threads: get_resources("fastqc_raw")["threads"]
    resources:
        mem_mb=get_resources("fastqc_raw")["mem_mb"],
        runtime=get_resources("fastqc_raw")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/fastqc_raw/{sample}_R{read}.txt",
    log:
        "results/pipeline_info/logs/fastqc/raw/{sample}_R{read}.log",
    conda:
        FASTQC_ENV
    shell:
        # shell:, not run:. Snakemake does NOT activate a rule's conda env for
        # shell() calls inside a run: block, so the conda: above was silently
        # doing nothing -- this only ever worked because fastqc also happens to
        # be in the monolithic workflow environment, and it would have broken
        # under `--sdm conda`. A run: block also pins the job to the main
        # snakemake process, so it could not be dispatched as a cluster job.
        #
        # The output .zip is named directly by the symlink's basename, so the
        # rename this rule used to do afterwards was provably a no-op and is
        # gone. trap EXIT removes the symlink even when FastQC fails.
        "mkdir -p {params.outdir} && "
        "rm -f {params.canonical} && "
        "ln -s {params.fq_abs} {params.canonical} && "
        "trap 'rm -f {params.canonical}' EXIT && "
        "fastqc -t {threads} -o {params.outdir} {params.canonical} > {log} 2>&1"


rule benchmark_summary:
    # Aggregates every benchmark file this configuration produces into a
    # self-describing MultiQC custom content table ("Resource usage",
    # id resource_usage), rendered after the RSeQC section inside the
    # multiqc_report.html via module_order in multiqc_config.yaml.
    input:
        all_benchmark_files(),
    output:
        "results/pipeline_info/benchmark_summary_mqc.json",
    params:
        allocated=allocated_resources_by_rule(),
    threads: get_resources("benchmark_summary")["threads"]
    resources:
        mem_mb=get_resources("benchmark_summary")["mem_mb"],
        runtime=get_resources("benchmark_summary")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/benchmark_summary/benchmark_summary.txt",
    log:
        "results/pipeline_info/logs/multiqc/benchmark_summary.log",
    conda:
        MULTIQC_ENV
    script:
        "../scripts/benchmark_summary.py"


def _chimera_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the chimera sample-QC view (rendered
    as interactive PCA + sample-distance plots inside the report). Only when
    the chimera stage is enabled and a counts matrix is written. Their
    directory is added to the MultiQC scan dirs via params.indirs."""
    if not CHIMERA_JUNCTION_ENABLED:
        return []
    if not config["chimera"]["junction"]["outputs"]["write_counts_matrix"]:
        return []
    transform = config["chimera"]["junction"]["qc"]["pca_transform"]
    return [
        f"results/chimera/qc/pca_{transform}_mqc.json",
        f"results/chimera/qc/heatmap_{transform}_mqc.json",
    ]


def _junction_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the chimera junction-QC barplot
    (per-sample direction composition) and the gene-TE chimeras barplot (the
    gene<->TE subset). Only when the chimera stage is enabled; independent of
    write_counts_matrix, since junction QC runs for every sample the annotator
    produces."""
    if not CHIMERA_JUNCTION_ENABLED:
        return []
    return [
        "results/chimera/qc/junction_qc_mqc.json",
        "results/chimera/qc/te_gene_chimeras_mqc.json",
        "results/chimera/qc/canonical_rate_mqc.json",
        "results/chimera/qc/junction_highlights_mqc.json",
        "results/chimera/qc/chimera_evidence_correlation_mqc.json",
        "results/chimera/qc/chimera_evidence_candidates_mqc.json",
        "results/chimera/qc/sample_evidence_status_mqc.json",
    ]


def _tecount_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the TEcounts sample-QC view (rendered
    as interactive PCA + sample-distance plots inside the report). Only when
    tetranscripts.qc.enabled (the default)."""
    if not TECOUNT_QC_ENABLED:
        return []
    transform = TECOUNT_QC["pca_transform"]
    return [
        f"results/tecount/qc/pca_{transform}_mqc.json",
        f"results/tecount/qc/heatmap_{transform}_mqc.json",
    ]


def _tecount_summary_mqc_inputs():
    """MultiQC custom-content JSONs for the per-sample summary barplots
    (gene-vs-TE assignment + TE class composition). Always produced, since
    they only need the raw cntTables -- independent of
    tetranscripts.qc.enabled."""
    return [
        "results/tecount/qc/tecount_assignment_mqc.json",
        "results/tecount/qc/tecount_te_class_mqc.json",
    ]


def _chimera_assembly_mqc_inputs():
    """MultiQC custom-content JSONs from the chimera-assembly screen (candidates
    by class + the "what to look at" highlights). Only when chimera.assembly
    is enabled (off by default)."""
    if not CHIMERA_ASSEMBLY_ENABLED:
        return []
    return [
        "results/chimera/qc/chimera_assembly_classes_mqc.json",
        "results/chimera/qc/chimera_assembly_highlights_mqc.json",
    ]


def _telocal_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the TElocal sample-QC view (rendered
    as interactive PCA + sample-distance plots inside the report) plus the
    per-sample summary barplots. The QC view only when telocal.qc.enabled
    (the default); the barplots always, since they only need the raw
    cntTables."""
    if not TELOCAL_ENABLED:
        return []
    files = [
        "results/telocal/qc/telocal_assignment_mqc.json",
        "results/telocal/qc/telocal_te_class_mqc.json",
    ]
    if TELOCAL_QC_ENABLED:
        transform = TELOCAL_QC["pca_transform"]
        files += [
            f"results/telocal/qc/pca_{transform}_mqc.json",
            f"results/telocal/qc/heatmap_{transform}_mqc.json",
        ]
    return files


rule evidence_overview:
    # The report's "TE analysis" section: which evidence layers this run
    # actually produced and which of them are independent of each other.
    # Reads the resolved switches from params (like config_used) rather than
    # re-deriving them, so it cannot drift from what really ran.
    output:
        "results/pipeline_info/evidence_overview_mqc.json",
    threads: get_resources("evidence_overview")["threads"]
    resources:
        mem_mb=get_resources("evidence_overview")["mem_mb"],
        runtime=get_resources("evidence_overview")["runtime"],
    log:
        "results/pipeline_info/logs/multiqc/evidence_overview.log",
    params:
        _sample_count=len(SAMPLES),
        _has_condition=HAS_CONDITION,
        _telocal_enabled=TELOCAL_ENABLED,
        _chimera_junction_enabled=CHIMERA_JUNCTION_ENABLED,
        _chimera_assembly_enabled=CHIMERA_ASSEMBLY_ENABLED,
        _two_pass=STAR_TWO_PASS,
    script:
        "../scripts/evidence_overview_mqc.py"


rule config_used:
    # The resolved run config, written as a MultiQC custom-content table so
    # the report records exactly which settings were used.
    output:
        "results/pipeline_info/config_used_mqc.json",
    threads: get_resources("config_used")["threads"]
    resources:
        mem_mb=get_resources("config_used")["mem_mb"],
        runtime=get_resources("config_used")["runtime"],
    log:
        "results/pipeline_info/logs/multiqc/config_used.log",
    params:
        _samples=SAMPLES,
        _sample_count=len(SAMPLES),
        _sjdb_overhang=SJDB_OVERHANG,
        _star_index=STAR_INDEX_DIR,
        _trim_enabled=TRIM_ENABLED,
        _tecount_qc_enabled=TECOUNT_QC_ENABLED,
        _tecount_qc=TECOUNT_QC,
        _chimera_enabled=CHIMERA_JUNCTION_ENABLED,
        _telocal_enabled=TELOCAL_ENABLED,
        _telocal_locind_auto=not _telocal_locind_cfg,
        _telocal_qc_enabled=TELOCAL_QC_ENABLED,
        _telocal_qc=TELOCAL_QC,
        _keep_merged_fastq=KEEP_MERGED_FASTQ,
        _keep_trimmed_fastq=KEEP_TRIMMED_FASTQ,
        _keep_star_index=KEEP_STAR_INDEX,
        _keep_telocal_index=KEEP_TELOCAL_INDEX,
    script:
        "../scripts/config_used_mqc.py"


rule multiqc:
    # Aggregates STAR alignment logs, (if trimming is enabled) TrimGalore!
    # + FastQC reports, the always-on raw FastQC reports, (if strandedness
    # auto-detection was used) RSeQC infer_experiment.py reports, the
    # always-on samtools flagstat + RSeQC read_distribution/geneBody_coverage
    # reports (all samples, not just AUTO_SAMPLES), the chimera
    # sample-QC (PCA + sample distances) and junction-QC barplot, the TEcounts
    # sample-QC (gated on tetranscripts.qc.enabled) and the always-on
    # per-sample summary barplots, the TElocal summary barplots + sample-QC
    # view (the latter gated on telocal.qc.enabled), the per-rule
    # benchmark/resource summary, and the pinned tool versions into one HTML
    # report. Runs the MultiQC version
    # pinned in config["versions"]["multiqc"]. The custom config
    # (multiqc_config.yaml) trims "_val_1"/"_val_2"/"_trimmed" off the FastQC
    # sample names so every module's rows merge into one clean row per sample
    # in General Stats.
    input:
        expand("results/star/{sample}_Log.final.out", sample=SAMPLES),
        expand("results/rseqc/{sample}_infer_experiment.txt",
               sample=STRAND_CHECK_SAMPLES),
        "results/rseqc/strandedness_check_mqc.json",
        expand("results/samtools/{sample}_flagstat.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}_read_distribution.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        all_fastqc_reports(),
        all_raw_fastqc_reports(),
        "results/pipeline_info/benchmark_summary_mqc.json",
        "results/pipeline_info/config_used_mqc.json",
        "results/pipeline_info/evidence_overview_mqc.json",
        "results/versions/rnaseq_mqc_versions.yml",
        chimera_qc=_chimera_qc_mqc_inputs(),
        junction_qc=_junction_qc_mqc_inputs(),
        tecount_qc=_tecount_qc_mqc_inputs(),
        tecount_summary=_tecount_summary_mqc_inputs(),
        telocal=_telocal_qc_mqc_inputs(),
        chimera_assembly=_chimera_assembly_mqc_inputs(),
        # Tracked so that editing it rebuilds the report (see common.smk).
        multiqc_config=MULTIQC_CONFIG,
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_report_data"),
    params:
        extra="",
        # Search directories come from the data inputs only. The config is an
        # input so its edits are tracked, but its directory holds the
        # workflow's own default-config YAMLs and must never be scanned.
        indirs=lambda wc, input: sorted(
            {os.path.dirname(f) for f in input if f != MULTIQC_CONFIG}
        ),
    threads: get_resources("multiqc")["threads"]
    resources:
        mem_mb=get_resources("multiqc")["mem_mb"],
        runtime=get_resources("multiqc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/multiqc/multiqc.txt",
    log:
        "results/pipeline_info/logs/multiqc/multiqc.log",
    conda:
        MULTIQC_ENV
    shell:
        # On failure, echo the log tail to stderr. MultiQC writes everything
        # to {log}, so a CI failure otherwise shows only "check log file(s)
        # for error details" -- and on a hosted runner that file is gone
        # with the workspace, which made several real failures undiagnosable
        # from the run output alone.
        "( multiqc {params.extra} --force "
        "-c {input.multiqc_config} "
        "-o results/qc -n multiqc_report "
        "{params.indirs} "
        "> {log} 2>&1 ) || "
        "( echo '-- multiqc failed; log tail --' >&2; "
        "tail -n 60 {log} >&2; exit 1 )\n"
