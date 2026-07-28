rule star_index:
    # Builds the STAR genome index once; reused by every sample's alignment.
    # Uses the STAR version pinned in config["versions"]["star"] (see the
    # generated env in common.smk) rather than a snakemake-wrapper, so the
    # exact tool version is fully under your control.
    input:
        fasta=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"],
    output:
        directory(config["star"]["index"]),
    params:
        sjdb_overhang=config["ref"]["sjdb_overhang"],
        extra="",
    threads: 8
    log:
        "logs/star/index.log",
    conda:
        STAR_ENV
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--runThreadN {threads} "
        "{params.extra} "
        "> {log} 2>&1"


rule gtf_to_genepred:
    # RSeQC's infer_experiment.py needs a BED12 reference model, not a GTF.
    # UCSC's gtfToGenePred + genePredToBed is the standard conversion route.
    input:
        gtf=config["ref"]["gtf"],
    output:
        genepred=temp("resources/annotation.genePred"),
    log:
        "logs/rseqc/gtf_to_genepred.log",
    conda:
        UCSC_TOOLS_ENV
    shell:
        "gtfToGenePred -genePredExt -ignoreGroupsWithoutExons "
        "{input.gtf} {output.genepred} > {log} 2>&1"


rule genepred_to_bed12:
    input:
        genepred="resources/annotation.genePred",
    output:
        bed12="resources/annotation.bed12",
    log:
        "logs/rseqc/genepred_to_bed12.log",
    conda:
        UCSC_TOOLS_ENV
    shell:
        "genePredToBed {input.genepred} {output.bed12} > {log} 2>&1"
