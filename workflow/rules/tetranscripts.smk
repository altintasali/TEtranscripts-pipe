rule tecount:
    # Per-sample gene + TE quantification. No official snakemake-wrapper
    # exists for TEtranscripts/TEcount, so this runs the tool directly in a
    # conda env generated from config["versions"] (see common.smk).
    input:
        bam=quant_bam_input,
        gtf=GTF,
        te_gtf=TE_GTF,
        strandedness=strandedness_input,
    output:
        "results/tecount/{sample}.cntTable.gz",
    params:
        stranded=get_strandedness_param,
        mode=config["tetranscripts"]["mode"],
        extra=config["tetranscripts"]["extra"],
        outdir="results/tecount",
    threads: get_resources("tecount")["threads"]
    resources:
        mem_mb=get_resources("tecount")["mem_mb"],
        runtime=get_resources("tecount")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tecount/{sample}.txt",
    log:
        "results/pipeline_info/logs/tecount/{sample}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        # Whole chain wrapped in ( ... ) so the redirect captures TEcount's
        # own stdout/stderr too, not just gzip's -- `cmd1 && cmd2 > log`
        # only redirects cmd2, silently dropping TEcount's own messages.
        "(mkdir -p {params.outdir} && "
        "TEcount --format BAM --mode {params.mode} "
        "-b {input.bam} "
        "--GTF {input.gtf} --TE {input.te_gtf} "
        "--stranded {params.stranded} "
        "--project {wildcards.sample} "
        "--outdir {params.outdir} "
        "{params.extra} "
        "&& gzip -f {params.outdir}/{wildcards.sample}.cntTable) "
        "> {log} 2>&1"
