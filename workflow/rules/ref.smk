rule gunzip_reference:
    # Decompresses a gzipped reference file (fasta/gtf/te_gtf) so downstream
    # tools that expect plain text (STAR's --genomeFastaFiles/--sjdbGTFfile,
    # TEcount/TEtranscripts' --GTF/--TE) can consume it directly. Only
    # triggered for whichever of fasta/gtf/te_gtf are actually given as
    # .gz in config.yaml -- see REFERENCE_GZ_SOURCES in common.smk.
    input:
        lambda wc: REFERENCE_GZ_SOURCES[wc.stem],
    output:
        "resources/decompressed/{stem}",
    threads: get_resources("gunzip_reference")["threads"]
    resources:
        mem_mb=get_resources("gunzip_reference")["mem_mb"],
        runtime=get_resources("gunzip_reference")["runtime"],
    log:
        "logs/gunzip/{stem}.log",
    shell:
        "gunzip -c {input} > {output} 2> {log}"


rule star_index:
    # Builds the STAR genome index once; reused by every sample's alignment.
    # sjdbOverhang defaults to being auto-detected from the sample sheet's
    # fastq read lengths (see SJDB_OVERHANG in common.smk) unless
    # ref.sjdb_overhang is set to an explicit integer in config.yaml.
    # Uses the STAR version pinned in config["versions"]["star"] (see the
    # generated env in common.smk) rather than a snakemake-wrapper, so the
    # exact tool version is fully under your control.
    input:
        fasta=FASTA,
        gtf=GTF,
    output:
        directory(config["star"]["index"]),
    params:
        sjdb_overhang=SJDB_OVERHANG,
        extra=config["star"].get("index_extra", ""),
    threads: get_resources("star_index")["threads"]
    resources:
        mem_mb=get_resources("star_index")["mem_mb"],
        runtime=get_resources("star_index")["runtime"],
    log:
        "logs/star/index.log",
    conda:
        STAR_ENV
    shell:
        # STAR has a long-standing, still-unresolved crash-on-exit bug
        # ("double free or corruption" / segfault while freeing memory
        # *after* all output has already been written and closed) --
        # confirmed benign by STAR's own author:
        # https://groups.google.com/g/rna-star/c/3_ckDieghws
        # ("this problem is happening after STAR finished all
        # calculations, so it does not affect the results"). So if STAR
        # exits non-zero, don't trust that alone -- check whether the core
        # index files actually landed on disk before treating it as a
        # real failure.
        "mkdir -p {output} && "
        "(STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--runThreadN {threads} "
        "{params.extra} "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether the index was "
        "actually written successfully anyway (see rule comment re: known "
        "benign STAR exit-time crash)' >> {log}; "
        "test -s {output}/SA && test -s {output}/SAindex && "
        "test -s {output}/Genome))"


rule gtf_to_genepred:
    # RSeQC's infer_experiment.py needs a BED12 reference model, not a GTF.
    # UCSC's gtfToGenePred + genePredToBed is the standard conversion route.
    input:
        gtf=GTF,
    output:
        genepred=temp("resources/annotation.genePred"),
    threads: get_resources("gtf_to_genepred")["threads"]
    resources:
        mem_mb=get_resources("gtf_to_genepred")["mem_mb"],
        runtime=get_resources("gtf_to_genepred")["runtime"],
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
    threads: get_resources("genepred_to_bed12")["threads"]
    resources:
        mem_mb=get_resources("genepred_to_bed12")["mem_mb"],
        runtime=get_resources("genepred_to_bed12")["runtime"],
    log:
        "logs/rseqc/genepred_to_bed12.log",
    conda:
        UCSC_TOOLS_ENV
    shell:
        "genePredToBed {input.genepred} {output.bed12} > {log} 2>&1"
