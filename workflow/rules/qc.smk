rule multiqc:
    # Aggregates STAR alignment logs, (if trimming is enabled) TrimGalore!
    # + FastQC reports, and (if strandedness auto-detection was used) RSeQC
    # infer_experiment.py reports into one HTML report. Runs the MultiQC
    # version pinned in config["versions"]["multiqc"].
    input:
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES),
        expand("results/rseqc/{sample}/infer_experiment.txt", sample=AUTO_SAMPLES),
        all_trim_outputs(),
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
        "-o results/qc -n multiqc_report "
        "{params.indirs} "
        "> {log} 2>&1\n"