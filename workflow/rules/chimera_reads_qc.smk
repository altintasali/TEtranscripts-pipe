# -----------------------------------------------------------------------------
# Chimera sample-QC: PCA / sample-distance views of the chimera counts
# matrix (results/chimera/counts_matrix.tsv), DESeq2-normalized
# (vst/rlog) or log2. The views are written as MultiQC custom-content JSON
# (pca_{transform}_mqc.json, heatmap_{transform}_mqc.json) and rendered
# interactively inside multiqc_report.html.
#
# Rules:
#   chimera_reads_sample_qc_transform  normalize the counts matrix for the QC view
#   chimera_reads_sample_qc            PCA + sample-distance plots from the transformed matrix
#
# QC filters (chimera.qc: min_samples_present / min_total_counts / min_events
# / pca_transform) apply ONLY to this view -- the all-events catalog and the
# counts matrix are never reduced. Both rules run in the dedicated R env
# (TETRANSCRIPTS_ENV, common.smk): deseq2.
# -----------------------------------------------------------------------------
rule chimera_reads_sample_qc_transform:
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/sample_qc.R",
        counts="results/chimera/counts_matrix.tsv.gz",
    output:
        "results/chimera/qc/{transform}_counts.tsv.gz",
    params:
        samples=config["samples"],
        min_samples_present=CHIMERA_QC["min_samples_present"],
        min_total_counts=CHIMERA_QC["min_total_counts"],
    threads: get_resources("chimera_reads_sample_qc_transform")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_reads_sample_qc_transform"),
        runtime=get_resources("chimera_reads_sample_qc_transform")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_sample_qc_transform/{transform}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_reads_sample_qc/transform_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {input.script} "
        "--transform chimera {input.counts} {params.samples} {wildcards.transform} "
        "{params.min_samples_present} {params.min_total_counts} "
        "{output} > {log} 2>&1"


rule chimera_reads_sample_qc:
    # PCA scatter + sample-to-sample distance heatmap of the transformed
    # counts, colored by condition (sample sheet's "condition" column; absent
    # -> one "all" group), emitted as MultiQC custom-content JSON.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/sample_qc.R",
        transformed="results/chimera/qc/{transform}_counts.tsv.gz",
    output:
        pca="results/chimera/qc/pca_{transform}_mqc.json",
        heatmap="results/chimera/qc/heatmap_{transform}_mqc.json",
    params:
        samples=config["samples"],
        min_events=CHIMERA_QC["min_events"],
    threads: get_resources("chimera_reads_sample_qc")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_reads_sample_qc"),
        runtime=get_resources("chimera_reads_sample_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_sample_qc/{transform}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_reads_sample_qc/plots_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {input.script} "
        "--plots chimera {input.transformed} {params.samples} {params.min_events} "
        "{wildcards.transform} {output.pca} {output.heatmap} > {log} 2>&1"
