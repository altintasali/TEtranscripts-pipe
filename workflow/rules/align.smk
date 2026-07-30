rule star_align:
    # Aligns one sample against the shared genome index. Output BAM is left
    # "Unsorted" (STAR's native read order) on purpose: TEtranscripts/TEcount
    # require either unsorted or queryname-sorted input, see
    # https://github.com/mhammell-laboratory/TEtranscripts#recommendations-for-tetranscripts-input-files
    # Runs the STAR version pinned in config["versions"]["star"].
    # Restarted once on failure (common transient issue on shared clusters).
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
    threads: get_resources("star_align")["threads"]
    resources:
        mem_mb=get_resources("star_align")["mem_mb"],
        runtime=get_resources("star_align")["runtime"],
        restart=2,
    log:
        "logs/star/align/{sample}.log",
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
    # The -m value is derived from the rule's mem_mb divided across threads.
    input:
        "results/star/{sample}/Aligned.out.bam",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        extra=lambda wc: f"-m {get_resources('samtools_sort')['mem_mb'] // get_resources('samtools_sort')['threads']}M",
    threads: get_resources("samtools_sort")["threads"]
    resources:
        mem_mb=get_resources("samtools_sort")["mem_mb"],
        runtime=get_resources("samtools_sort")["runtime"],
    log:
        "logs/samtools/sort/{sample}.log",
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
    threads: get_resources("samtools_index")["threads"]
    resources:
        mem_mb=get_resources("samtools_index")["mem_mb"],
        runtime=get_resources("samtools_index")["runtime"],
    log:
        "logs/samtools/index/{sample}.log",
    shell:
        "samtools index {params.extra} -@ {threads} {input} {output} > {log} 2>&1"
