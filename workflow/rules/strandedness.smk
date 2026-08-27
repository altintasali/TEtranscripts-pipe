rule rseqc_infer_experiment:
    # Infers library strandedness from a sample of aligned reads against the
    # BED12 gene model. Runs the RSeQC version pinned in
    # config["versions"]["rseqc"].
    #
    # Requested for every sample in STRAND_CHECK_SAMPLES (common.smk): the
    # "auto" ones, which need the inference, plus -- under
    # strandedness.check_provided, the default -- the ones that declared a
    # value, so the report can flag a sample sheet that disagrees with the
    # data. For a declared sample the result is REPORTED ONLY; the sample
    # sheet still decides what TEcount is run with.
    input:
        aln="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        refgene="results/reference/annotation.bed12",
    output:
        "results/rseqc/{sample}_infer_experiment.txt",
    params:
        extra="",
    threads: get_resources("rseqc_infer_experiment")["threads"]
    resources:
        mem_mb=get_resources("rseqc_infer_experiment")["mem_mb"],
        runtime=get_resources("rseqc_infer_experiment")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/rseqc_infer_experiment/{sample}.txt",
    log:
        "results/pipeline_info/logs/rseqc/infer_experiment/{sample}.log",
    conda:
        RSEQC_ENV
    shell:
        "infer_experiment.py {params.extra} "
        "--input-file {input.aln} "
        "--refgene {input.refgene} "
        "> {output} 2> {log}"


rule determine_strandedness:
    # Converts the RSeQC infer_experiment.py report into the TEtranscripts/
    # TEcount --stranded value (no/forward/reverse) for this sample.
    input:
        txt="results/rseqc/{sample}_infer_experiment.txt",
    output:
        txt="results/rseqc/{sample}_strandedness.txt",
    params:
        min_fraction=config["strandedness"]["min_fraction"],
    threads: get_resources("determine_strandedness")["threads"]
    resources:
        mem_mb=get_resources("determine_strandedness")["mem_mb"],
        runtime=get_resources("determine_strandedness")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/determine_strandedness/{sample}.txt",
    log:
        "results/pipeline_info/logs/rseqc/determine_strandedness/{sample}.log",
    script:
        "../scripts/determine_strandedness.py"


rule strandedness_check:
    # Compares each sample's declared strandedness against what RSeQC
    # inferred, as a MultiQC table with mismatches flagged
    # (strandedness_check_mqc.py). Modelled on nf-core/rnaseq's strandedness
    # check.
    #
    # Worth a dedicated section because strandedness is the one setting that
    # is both easy to get wrong (kit docs say "stranded" without saying which
    # direction) and completely silent when wrong -- no job fails, the counts
    # are just wrong.
    input:
        reports=expand(
            "results/rseqc/{sample}_infer_experiment.txt",
            sample=STRAND_CHECK_SAMPLES,
        ),
        calls=expand(
            "results/rseqc/{sample}_strandedness.txt",
            sample=STRAND_CHECK_SAMPLES,
        ),
    output:
        "results/rseqc/strandedness_check_mqc.json",
    params:
        samples=STRAND_CHECK_SAMPLES,
        declared={s: SAMPLE_STRANDED_MODE[s] for s in STRAND_CHECK_SAMPLES},
    threads: get_resources("strandedness_check")["threads"]
    resources:
        mem_mb=get_resources("strandedness_check")["mem_mb"],
        runtime=get_resources("strandedness_check")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/strandedness_check/strandedness_check.txt",
    log:
        "results/pipeline_info/logs/rseqc/strandedness_check.log",
    script:
        "../scripts/strandedness_check_mqc.py"
