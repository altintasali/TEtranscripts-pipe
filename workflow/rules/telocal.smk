rule telocal_locind:
    # Build the TElocal locus index from the TE GTF.  Only triggered when
    # the user does not provide a pre-built .locInd via config
    # telocal.locind -- the telocal rule selects the auto-built path when
    # locind is empty.
    input:
        te_gtf=TE_GTF,
    output:
        "results/telocal.locInd",
    threads: get_resources("telocal_locind")["threads"]
    resources:
        mem_mb=get_resources("telocal_locind")["mem_mb"],
        runtime=get_resources("telocal_locind")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_locind/locind.txt",
    log:
        "results/pipeline_info/logs/telocal/locind.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "python3 {SCRIPTS_DIR}/build_telocal_index.py "
        "--gtf {input.te_gtf} "
        "--out {output} > {log} 2>&1"


rule telocal:
    # Per-sample locus-level TE quantification (TElocal). Complements TEcount's
    # subfamily-level quantification by resolving TEs per genomic instance.
    # Uses the same unsorted BAM as TEcount; no re-alignment needed.
    # TElocal has no --outdir flag, so we cd into the output directory before
    # running it; {input.bam} resolves to an absolute path so the cd is safe.
    input:
        bam="results/star/{sample}_Aligned.out.bam",
        gtf=GTF,
        locind=_telocal_locind_path(),
        strandedness=strandedness_input,
    output:
        "results/telocal/{sample}.cntTable.gz",
    params:
        stranded=get_strandedness_param,
        extra=config["telocal"]["extra"],
    threads: get_resources("telocal")["threads"]
    resources:
        mem_mb=get_resources("telocal")["mem_mb"],
        runtime=get_resources("telocal")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal/{sample}.txt",
    log:
        "results/pipeline_info/logs/telocal/{sample}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "( mkdir -p results/telocal &&"
        " TElocal --sortByPos -b {input.bam}"
        " --GTF {input.gtf} --TE {input.locind}"
        " --stranded {params.stranded}"
        " --project {wildcards.sample}"
        " {params.extra}"
        " && mv {wildcards.sample}.cntTable results/telocal/"
        " && gzip -f results/telocal/{wildcards.sample}.cntTable )"
        " > {log} 2>&1"


def telocal_counts_input():
    return [
        f"results/telocal/{s}.cntTable.gz"
        for s in SAMPLES
    ]


rule telocal_summary:
    # Per-sample TElocal summary stats for the MultiQC report (custom
    # content): gene-vs-TE assignment and TE class composition as counts and
    # percentages (telocal_summary_mqc.py). Uses the RAW cntTables.
    input:
        tables=telocal_counts_input(),
    output:
        assignment="results/telocal/qc/telocal_assignment_mqc.json",
        te_class="results/telocal/qc/telocal_te_class_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("telocal")["threads"]
    resources:
        mem_mb=get_resources("telocal")["mem_mb"],
        runtime=get_resources("telocal")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_summary/telocal_summary.txt",
    log:
        "results/pipeline_info/logs/telocal/summary.log",
    shell:
        "python3 {SCRIPTS_DIR}/telocal_summary_mqc.py "
        "--tables {input.tables} "
        "--samples {params.samples} "
        "--out-assignment {output.assignment} "
        "--out-class {output.te_class} > {log} 2>&1"
