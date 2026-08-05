rule tecount:
    # Per-sample gene + TE quantification. No official snakemake-wrapper
    # exists for TEtranscripts/TEcount, so this runs the tool directly in a
    # conda env generated from config["versions"] (see common.smk).
    input:
        bam="results/star/{sample}_Aligned.out.bam",
        gtf=GTF,
        te_gtf=TE_GTF,
        strandedness=strandedness_input,
    output:
        "results/tecount/{sample}.cntTable",
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
        "mkdir -p {params.outdir} && "
        "TEcount --format BAM --mode {params.mode} "
        "-b {input.bam} "
        "--GTF {input.gtf} --TE {input.te_gtf} "
        "--stranded {params.stranded} "
        "--project {wildcards.sample} "
        "--outdir {params.outdir} "
        "{params.extra} "
        "> {log} 2>&1"


rule tetranscripts_diffexp:
    # Differential gene/TE enrichment between treatment and control BAM
    # groups. Contrasts are derived automatically (see CONTRASTS in
    # common.smk) from every pairwise combination of distinct values in the
    # sample sheet's "condition" column; if that column is absent, CONTRASTS
    # is empty and this rule is simply never requested. Internally
    # re-quantifies all samples together and runs DESeq2.
    input:
        treatment=lambda wc: expand(
            "results/star/{sample}_Aligned.out.bam",
            sample=CONTRASTS[wc.contrast]["treatment"],
        ),
        control=lambda wc: expand(
            "results/star/{sample}_Aligned.out.bam",
            sample=CONTRASTS[wc.contrast]["control"],
        ),
        strandedness=contrast_strandedness_input,
        gtf=GTF,
        te_gtf=TE_GTF,
    output:
        cnt_table="results/tetranscripts/{contrast}.cntTable",
        deseq_script="results/tetranscripts/{contrast}_DESeq2.R",
        full="results/tetranscripts/{contrast}_gene_TE_analysis.txt",
        sig="results/tetranscripts/{contrast}_sigdiff_gene_TE.txt",
    params:
        stranded=get_contrast_strandedness_param,
        mode=config["tetranscripts"]["mode"],
        padj=config["tetranscripts"]["padj"],
        foldchange=config["tetranscripts"]["foldchange"],
        minread=config["tetranscripts"]["minread"],
        extra=config["tetranscripts"]["extra"],
        outdir="results/tetranscripts",
    threads: get_resources("tetranscripts_diffexp")["threads"]
    resources:
        mem_mb=get_resources("tetranscripts_diffexp")["mem_mb"],
        runtime=get_resources("tetranscripts_diffexp")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/tetranscripts_diffexp/{contrast}.txt",
    log:
        "results/pipeline_info/logs/tetranscripts/{contrast}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "mkdir -p {params.outdir} && "
        "TEtranscripts --format BAM --mode {params.mode} "
        "-t {input.treatment} "
        "-c {input.control} "
        "--GTF {input.gtf} --TE {input.te_gtf} "
        "--stranded {params.stranded} "
        "--project {wildcards.contrast} "
        "--padj {params.padj} --foldchange {params.foldchange} "
        "--minread {params.minread} "
        "--outdir {params.outdir} "
        "{params.extra} "
        "> {log} 2>&1"
