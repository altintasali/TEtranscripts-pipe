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
    # input even when trimming is disabled. FastQC names its .zip after the
    # input file's basename, which only matches the canonical
    # {sample}_R{read} name for the merged results/fastq/ files -- a raw
    # single-lane path would produce a differently-named zip, so it is
    # renamed to the declared output here.
    input:
        fq=lambda wc: merged_fastq_path(wc.sample, int(wc.read)),
    output:
        zip="results/fastqc/raw/{sample}_R{read}_fastqc.zip",
    params:
        outdir="results/fastqc/raw",
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
    run:
        os.makedirs(params.outdir, exist_ok=True)
        shell(
            "fastqc -t {threads} -o {params.outdir} {input.fq} "
            "> {log} 2>&1"
        )
        produced = os.path.join(
            params.outdir, _fastqc_base_name(input.fq) + "_fastqc.zip"
        )
        if os.path.abspath(produced) != os.path.abspath(output.zip):
            os.replace(produced, output.zip)


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
    if not CHIMERA_ENABLED:
        return []
    if not config["chimera"]["outputs"]["write_counts_matrix"]:
        return []
    transform = config["chimera"]["qc"]["pca_transform"]
    return [
        f"results/chimera/qc/pca_{transform}_mqc.json.gz",
        f"results/chimera/qc/heatmap_{transform}_mqc.json.gz",
    ]


def _junction_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the chimera junction-QC barplot
    (per-sample direction composition) and the TE-chimeras barplot (the
    gene<->TE subset). Only when the chimera stage is enabled; independent of
    write_counts_matrix, since junction QC runs for every sample the annotator
    produces."""
    if not CHIMERA_ENABLED:
        return []
    return [
        "results/chimera/qc/junction_qc_mqc.json.gz",
        "results/chimera/qc/te_chimeras_mqc.json.gz",
    ]


def _tecount_qc_mqc_inputs():
    """MultiQC custom-content JSONs from the TEcounts sample-QC view (rendered
    as interactive PCA + sample-distance plots inside the report). Only when
    tetranscripts.qc.enabled (the default)."""
    if not TECOUNT_QC_ENABLED:
        return []
    transform = TECOUNT_QC["pca_transform"]
    return [
        f"results/tecount/qc/pca_{transform}_mqc.json.gz",
        f"results/tecount/qc/heatmap_{transform}_mqc.json.gz",
    ]


def _tecount_summary_mqc_inputs():
    """MultiQC custom-content JSONs for the per-sample summary barplots
    (gene-vs-TE assignment + TE class composition). Always produced, since
    they only need the raw cntTables -- independent of
    tetranscripts.qc.enabled."""
    return [
        "results/tecount/qc/tecount_assignment_mqc.json.gz",
        "results/tecount/qc/tecount_te_class_mqc.json.gz",
    ]


rule config_used:
    # The resolved run config, written as a MultiQC custom-content table so
    # the report records exactly which settings were used.
    output:
        "results/pipeline_info/config_used_mqc.json.gz",
    threads: get_resources("config_used")["threads"]
    resources:
        mem_mb=get_resources("config_used")["mem_mb"],
        runtime=get_resources("config_used")["runtime"],
    log:
        "results/pipeline_info/logs/multiqc/config_used.log",
    run:
        import gzip
        import json

        star_extra = config.get("star", {}).get("extra", "") or "(none)"
        trim_extra = config.get("trimming", {}).get("extra", "") or "(none)"
        te_extra = config.get("tetranscripts", {}).get("extra", "") or "(none)"

        rows = {
            "pipeline_version": open("VERSION").read().strip(),
            "samples": f"{len(SAMPLES)} ({', '.join(SAMPLES)})",
            "ref.fasta": str(config["ref"]["fasta"]),
            "ref.gtf": str(config["ref"]["gtf"]),
            "ref.te_gtf": str(config["ref"]["te_gtf"]),
            "ref.sjdb_overhang": str(SJDB_OVERHANG),
            "star.index": STAR_INDEX_DIR,
            "star.extra": star_extra,
            "trimming.enabled": str(TRIM_ENABLED),
            "trimming.trim_nextseq": str(config.get("trimming", {}).get("trim_nextseq", 0)),
            "trimming.extra": trim_extra,
            "strandedness.min_fraction": str(config.get("strandedness", {}).get("min_fraction", 0.6)),
            "tetranscripts.mode": config["tetranscripts"]["mode"],
            "tetranscripts.padj": str(config["tetranscripts"]["padj"]),
            "tetranscripts.foldchange": str(config["tetranscripts"]["foldchange"]),
            "tetranscripts.minread": str(config["tetranscripts"]["minread"]),
            "tetranscripts.extra": te_extra,
            "tetranscripts.qc.enabled": str(TECOUNT_QC_ENABLED),
            "tetranscripts.qc.feature_class": TECOUNT_QC["feature_class"],
            "tetranscripts.qc.pca_transform": TECOUNT_QC["pca_transform"],
            "chimera.enabled": str(CHIMERA_ENABLED),
            "chimera.breakpoint_tolerance": str(config.get("chimera", {}).get("breakpoint_tolerance", 0)),
            "chimera.require_canonical_junction": str(config.get("chimera", {}).get("require_canonical_junction", False)),
            "chimera.qc.pca_transform": config.get("chimera", {}).get("qc", {}).get("pca_transform", "vst"),
            "outputs.keep_merged_fastq": str(KEEP_MERGED_FASTQ),
            "outputs.keep_trimmed_fastq": str(KEEP_TRIMMED_FASTQ),
        }

        doc = {
            "id": "config_used",
            "section_name": "Configuration used",
            "description": (
                "The resolved run configuration (config values from the "
                "config file passed with --configfile)."
            ),
            "plot_type": "table",
            "pconfig": {
                "id": "config_used_table",
                "title": "Configuration used",
                "col1_header": "Setting",
                "sort_rows": False,
            },
            "headers": {"value": {"title": "Value"}},
            "data": {k: {"value": v} for k, v in rows.items()},
        }
        with gzip.open(str(output), "wt") as fh:
            json.dump(doc, fh, indent=2)


rule multiqc:
    # Aggregates STAR alignment logs, (if trimming is enabled) TrimGalore!
    # + FastQC reports, the always-on raw FastQC reports, (if strandedness
    # auto-detection was used) RSeQC infer_experiment.py reports, the chimera
    # sample-QC (PCA + sample distances) and junction-QC barplot, the TEcounts
    # sample-QC (gated on tetranscripts.qc.enabled) and the always-on
    # per-sample summary barplots, the
    # per-rule benchmark/resource summary, and the pinned tool versions into
    # one HTML report. Runs the MultiQC version
    # pinned in config["versions"]["multiqc"]. The custom config
    # (multiqc_config.yaml) trims "_val_1"/"_val_2"/"_trimmed" off the FastQC
    # sample names so every module's rows merge into one clean row per sample
    # in General Stats.
    input:
        expand("results/star/{sample}_Log.final.out", sample=SAMPLES),
        expand("results/rseqc/{sample}_infer_experiment.txt", sample=AUTO_SAMPLES),
        all_fastqc_reports(),
        all_raw_fastqc_reports(),
        "results/pipeline_info/benchmark_summary_mqc.json",
        "results/pipeline_info/config_used_mqc.json.gz",
        "results/versions/rnaseq_mqc_versions.yml.gz",
        chimera_qc=_chimera_qc_mqc_inputs(),
        junction_qc=_junction_qc_mqc_inputs(),
        tecount_qc=_tecount_qc_mqc_inputs(),
        tecount_summary=_tecount_summary_mqc_inputs(),
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_report_data"),
    params:
        extra="",
        indirs=lambda wc, input: sorted({os.path.dirname(f) for f in input}),
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
        "multiqc {params.extra} --force "
        "-c workflow/default-config/multiqc_config.yaml "
        "-o results/qc -n multiqc_report "
        "{params.indirs} "
        "> {log} 2>&1\n"
