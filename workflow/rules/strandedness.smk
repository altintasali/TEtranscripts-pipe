rule rseqc_infer_experiment:
    # Infers library strandedness from a sample of aligned reads against the
    # BED12 gene model. Only needed when strandedness.mode == "auto". Runs
    # the RSeQC version pinned in config["versions"]["rseqc"].
    input:
        aln="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
        refgene="resources/annotation.bed12",
    output:
        "results/rseqc/{sample}/infer_experiment.txt",
    params:
        extra="",
    log:
        "logs/rseqc/infer_experiment/{sample}.log",
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
        txt="results/rseqc/{sample}/infer_experiment.txt",
    output:
        txt="results/rseqc/{sample}/strandedness.txt",
    params:
        min_fraction=config["strandedness"]["min_fraction"],
    log:
        "logs/rseqc/determine_strandedness/{sample}.log",
    script:
        "../scripts/determine_strandedness.py"
