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
    threads: get_resources("star_align")["threads"]
    resources:
        mem_mb=get_resources("star_align")["mem_mb"],
        runtime=get_resources("star_align")["runtime"],
    log:
        "logs/star/align/{sample}.log",
    conda:
        STAR_ENV
    shell:
        # See the star_index rule's comment: STAR has a known benign
        # crash-on-exit bug (confirmed by its author) that can happen
        # after alignment has already finished and all output files were
        # written -- so a non-zero exit is checked against actual output
        # completeness rather than trusted outright.
        "mkdir -p {params.prefix} && "
        "(STAR --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {params.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM Unsorted "
        "{params.extra} "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether alignment output "
        "is actually complete anyway (see rule comment re: known benign "
        "STAR exit-time crash)' >> {log}; "
        "test -s {output.aln} && test -s {output.log_final} && "
        "grep -q 'ALL DONE!' {output.log_final}))"


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
    threads: get_resources("samtools_sort")["threads"]
    resources:
        mem_mb=get_resources("samtools_sort")["mem_mb"],
        runtime=get_resources("samtools_sort")["runtime"],
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
    threads: get_resources("samtools_index")["threads"]
    resources:
        mem_mb=get_resources("samtools_index")["mem_mb"],
        runtime=get_resources("samtools_index")["runtime"],
    log:
        "logs/samtools/index/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools index {params.extra} -@ {threads} {input} {output} > {log} 2>&1"
