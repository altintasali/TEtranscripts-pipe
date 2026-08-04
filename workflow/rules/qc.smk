rule multiqc:
    # Aggregates STAR alignment logs, (if trimming is enabled) TrimGalore!
    # + FastQC reports, (if strandedness auto-detection was used) RSeQC
    # infer_experiment.py reports, and the pinned tool versions into one HTML
    # report. Runs the MultiQC version pinned in
    # config["versions"]["multiqc"]. The custom config (multiqc_config.yaml)
    # trims "_val_1"/"_val_2"/"_trimmed" off the FastQC sample names so every
    # module's rows merge into one clean row per sample in General Stats.
    input:
        expand("results/star/{sample}_Log.final.out", sample=SAMPLES),
        expand("results/rseqc/{sample}_infer_experiment.txt", sample=AUTO_SAMPLES),
        all_fastqc_reports(),
        "results/versions/rnaseq_mqc_versions.yml",
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
        "results/pipeline_info/benchmarks/multiqc.txt",
    log:
        "logs/multiqc.log",
    conda:
        MULTIQC_ENV
    shell:
        "multiqc {params.extra} --force "
        "-c workflow/default-config/multiqc_config.yaml "
        "-o results/qc -n multiqc_report "
        "{params.indirs} "
        "> {log} 2>&1\n"
