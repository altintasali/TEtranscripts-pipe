# -----------------------------------------------------------------------------
# TEcounts sample-QC: PCA / sample-distance views of the per-sample TEcount
# count tables (results/tecount/{sample}.cntTable), DESeq2-normalized
# (vst/rlog) or log2. TEcount itself quantifies genes + TEs into one table;
# tecount_counts merges those into results/tecount/counts_matrix.tsv.gz
# unfiltered (genes + TEs, exactly what TEcount reports), and a second,
# QC-only step (tecount_qc_counts) filters a copy of that matrix down to TE
# subfamilies (the default) or genes via tecounts.qc.feature_class -- only
# the PCA/clustering view sees the filtered copy.
#
# Rules:
#   tecount_counts      merge per-sample cntTables -> counts matrix (genes + TEs)
#   tecount_qc_counts   filter a copy of that matrix to the QC view's feature class
#   tecount_qc_transform normalize the filtered matrix for the QC view
#   tecount_qc          PCA + sample-distance plots from the transformed matrix
#   tecount_summary     per-sample assignment + TE-class summary barplots for
#                       the MultiQC report (raw cntTables, pure python)
#
# The view reuses the shared sample_qc.R script (same as the chimera sample-QC
# stage, passed --view tecount) and runs in the TETRANSCRIPTS_ENV, which
# already ships DESeq2 + R (the R env the chimera stage uses is a separate
# build of the same tools; reusing the TEtranscripts env keeps the download
# count down).
#
# QC filters (tetranscripts.qc: min_samples_present / min_total_counts /
# min_events / pca_transform / feature_class) apply ONLY to this view -- the
# per-sample cntTables and the top-level counts_matrix.tsv.gz are never
# reduced.
# -----------------------------------------------------------------------------


def all_tecount_qc_outputs():
    """TEcounts sample-QC artifacts for the `all` target (Snakefile). Only
    produced when tetranscripts.qc.enabled is true."""
    transform = TECOUNT_QC["pca_transform"]
    return [
        "results/tecount/counts_matrix.tsv.gz",
        "results/tecount/qc/counts_matrix.tsv.gz",
        f"results/tecount/qc/{transform}_counts.tsv.gz",
        f"results/tecount/qc/pca_{transform}_mqc.json",
        f"results/tecount/qc/heatmap_{transform}_mqc.json",
    ]


def all_tecount_summary_outputs():
    """Per-sample TEcounts summary-barplot JSONs (gene-vs-TE assignment + TE
    class composition) for the `all` target. Always produced: they only need
    the raw cntTables, so they are independent of tetranscripts.qc.enabled."""
    return [
        "results/tecount/qc/tecount_assignment_mqc.json",
        "results/tecount/qc/tecount_te_class_mqc.json",
    ]


def tecount_counts_input():
    return [
        f"results/tecount/{s}.cntTable.gz"
        for s in SAMPLES
    ]


rule tecount_counts:
    # Merges every sample's TEcount count table into the feature x sample
    # counts matrix (tecount_counts.py), always keeping every feature --
    # genes and TEs -- exactly as TEcount reports them. This is the native
    # merged deliverable; the QC-view's feature-class filtering happens
    # downstream, in tecount_qc_counts, on a separate copy.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/tecount_counts.py",
        tables=tecount_counts_input(),
    output:
        counts="results/tecount/counts_matrix.tsv.gz",
    params:
        sample_names=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("tecount_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("tecount_counts"),
        runtime=get_resources("tecount_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount_counts/tecount_counts.txt",
    log:
        "results/pipeline_info/logs/tecount/counts.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--feature-class all "
        "--out-counts {output.counts} > {log} 2>&1"


rule tecount_qc_counts:
    # Filters a copy of the full counts matrix (genes + TEs) down to the
    # QC-view's feature class -- TE subfamilies (default), genes, or all --
    # via tecounts.qc.feature_class. results/tecount/counts_matrix.tsv.gz
    # itself is untouched; only this QC-scoped copy is filtered, and only
    # the PCA/clustering view (tecount_qc_transform) reads it.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/tecount_counts.py",
        matrix="results/tecount/counts_matrix.tsv.gz",
        gtf=GTF,
        te_gtf=TE_GTF,
    output:
        "results/tecount/qc/counts_matrix.tsv.gz",
    params:
        feature_class=TECOUNT_QC["feature_class"],
    threads: get_resources("tecount_qc_counts")["threads"]
    resources:
        mem_mb=get_resources("tecount_qc_counts")["mem_mb"],
        runtime=get_resources("tecount_qc_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount_qc_counts/tecount_qc_counts.txt",
    log:
        "results/pipeline_info/logs/tecount/qc/counts.log",
    shell:
        "python3 {input.script} "
        "--in-matrix {input.matrix} "
        "--key-style tecount "
        "--gtf {input.gtf} --te-gtf {input.te_gtf} "
        "--feature-class {params.feature_class} "
        "--out-counts {output} > {log} 2>&1"


rule tecount_qc_transform:
    # Normalizes the TEcounts counts matrix for the QC view (sample_qc.R
    # --transform tecount). vst/rlog use DESeq2's blind normalization; log2 is
    # log2(x + 1) without DESeq2. Filters in tetranscripts.qc apply only to
    # this view.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/sample_qc.R",
        counts="results/tecount/qc/counts_matrix.tsv.gz",
    output:
        "results/tecount/qc/{transform}_counts.tsv.gz",
    params:
        samples=config["samples"],
        min_samples_present=TECOUNT_QC["min_samples_present"],
        min_total_counts=TECOUNT_QC["min_total_counts"],
    threads: get_resources("tecount_qc_transform")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("tecount_qc_transform"),
        runtime=get_resources("tecount_qc_transform")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount_qc_transform/{transform}.txt",
    log:
        "results/pipeline_info/logs/tecount/qc/transform_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {input.script} "
        "--transform tecount {input.counts} {params.samples} {wildcards.transform} "
        "{params.min_samples_present} {params.min_total_counts} "
        "{output} > {log} 2>&1"


rule tecount_qc:
    # PCA scatter + sample-to-sample distance heatmap of the transformed
    # TEcounts counts, colored by condition (sample sheet's "condition"
    # column; absent -> one "all" group), emitted as MultiQC custom-content
    # JSON (ids tecount_chimera_reads_sample_qc_pca / tecount_chimera_reads_sample_qc_heatmap, ordered
    # inside the custom_content module by multiqc_config.yaml).
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/sample_qc.R",
        transformed="results/tecount/qc/{transform}_counts.tsv.gz",
    output:
        pca="results/tecount/qc/pca_{transform}_mqc.json",
        heatmap="results/tecount/qc/heatmap_{transform}_mqc.json",
    params:
        samples=config["samples"],
        min_events=TECOUNT_QC["min_events"],
    threads: get_resources("tecount_qc")["threads"]
    resources:
        mem_mb=get_resources("tecount_qc")["mem_mb"],
        runtime=get_resources("tecount_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount_qc/{transform}.txt",
    log:
        "results/pipeline_info/logs/tecount/qc/plots_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {input.script} "
        "--plots tecount {input.transformed} {params.samples} {params.min_events} "
        "{wildcards.transform} {output.pca} {output.heatmap} > {log} 2>&1"


rule tecount_summary:
    # Per-sample TEcounts summary stats for the MultiQC report (custom
    # content): gene-vs-TE assignment and TE class composition as counts and
    # percentages (tecount_summary_mqc.py). Unlike the sample-QC view it uses
    # the RAW cntTables (all features), so it is independent of
    # tetranscripts.qc.enabled and .feature_class, and needs no R env --
    # always produced.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/tecount_summary_mqc.py",
        tables=tecount_counts_input(),
    output:
        assignment="results/tecount/qc/tecount_assignment_mqc.json",
        te_class="results/tecount/qc/tecount_te_class_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("tecount_summary")["threads"]
    resources:
        mem_mb=get_resources("tecount_summary")["mem_mb"],
        runtime=get_resources("tecount_summary")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount_summary/tecount_summary.txt",
    log:
        "results/pipeline_info/logs/tecount/qc/summary.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--samples {params.samples} "
        "--out-assignment {output.assignment} "
        "--out-class {output.te_class} > {log} 2>&1"
