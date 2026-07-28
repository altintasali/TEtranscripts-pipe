rule star_align:
    # Aligns one sample against the shared genome index. Output BAM is left
    # "Unsorted" (STAR's native read order) on purpose: TEtranscripts/TEcount
    # require either unsorted or queryname-sorted input, see
    # https://github.com/mhammell-laboratory/TEtranscripts#recommendations-for-tetranscripts-input-files
    # Runs the STAR version pinned in config["versions"]["star"].
    input:
        unpack(star_input),
        idx=config["star"]["index"],
    output:
        aln="results/star/{sample}/Aligned.out.bam",
        log_final="results/star/{sample}/Log.final.out",
        sj="results/star/{sample}/SJ.out.tab",
    params:
        reads=star_reads_param,
        read_command=star_read_command_param,
        prefix=lambda wc: f"results/star/{wc.sample}/",
        extra=config["star"]["extra"],
    threads: 12
    log:
        "logs/star/align/{sample}.log",
    conda:
        STAR_ENV
    shell:
        "mkdir -p {params.prefix} && "
        "STAR --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {params.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM Unsorted "
        "{params.extra} "
        "> {log} 2>&1"


rule samtools_sort:
    # A coordinate-sorted + indexed copy is only needed for RSeQC/QC purposes;
    # TEcount/TEtranscripts always consume the unsorted STAR output above.
    # Runs the samtools version pinned in config["versions"]["samtools"].
    input:
        "results/star/{sample}/Aligned.out.bam",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        extra="-m 3G",
    threads: 8
    log:
        "logs/samtools/sort/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools sort {params.extra} -@ {threads} "
        "-o {output} {input} > {log} 2>&1"


rule samtools_index:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    params:
        extra="",
    threads: 4
    log:
        "logs/samtools/index/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools index {params.extra} -@ {threads} {input} {output} > {log} 2>&1"
