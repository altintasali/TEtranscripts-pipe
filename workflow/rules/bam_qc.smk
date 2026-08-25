rule samtools_flagstat:
    # Feeds MultiQC's built-in Samtools module -- mapped %, properly-paired
    # %, duplicate count, secondary/supplementary read counts. Complements
    # STAR's own alignment-rate stats with samtools' read-flag breakdown.
    # Runs for every sample (not just AUTO_SAMPLES), unlike
    # rseqc_infer_experiment. Runs the samtools version pinned in
    # config["versions"]["samtools"].
    input:
        bam="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/samtools/{sample}_flagstat.txt",
    threads: get_resources("samtools_flagstat")["threads"]
    resources:
        mem_mb=get_resources("samtools_flagstat")["mem_mb"],
        runtime=get_resources("samtools_flagstat")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/samtools_flagstat/{sample}.txt",
    log:
        "results/pipeline_info/logs/samtools/flagstat/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output} 2> {log}"


rule rseqc_read_distribution:
    # Read distribution across genomic features (CDS/UTR/intron/intergenic)
    # -- MultiQC's RSeQC module renders this as its own subsection alongside
    # Infer experiment, so no multiqc_config.yaml change is needed beyond
    # the file existing. Flags library-prep issues (e.g. high
    # intronic/intergenic fraction suggesting DNA contamination or
    # incomplete polyA selection).
    input:
        aln="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        refgene="results/reference/annotation.bed12",
    output:
        "results/rseqc/{sample}_read_distribution.txt",
    threads: get_resources("rseqc_read_distribution")["threads"]
    resources:
        mem_mb=get_resources("rseqc_read_distribution")["mem_mb"],
        runtime=get_resources("rseqc_read_distribution")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/rseqc_read_distribution/{sample}.txt",
    log:
        "results/pipeline_info/logs/rseqc/read_distribution/{sample}.log",
    conda:
        RSEQC_ENV
    shell:
        "read_distribution.py -i {input.aln} -r {input.refgene} "
        "> {output} 2> {log}"


rule rseqc_gene_body_coverage:
    # 5'->3' gene body coverage -- flags RNA degradation (3' bias) or
    # incomplete reverse-transcription (5' bias); MultiQC's RSeQC module
    # picks up geneBody_coverage.py's .txt table as its own subsection.
    # geneBody_coverage.py names outputs from -o's prefix (.txt/.r/.pdf);
    # only the .txt table is a tracked output -- the .r/.pdf are RSeQC's own
    # plotting scaffolding, unused since MultiQC parses the .txt directly.
    input:
        aln="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        refgene="results/reference/annotation.bed12",
    output:
        "results/rseqc/{sample}.geneBodyCoverage.txt",
    params:
        prefix="results/rseqc/{sample}",
    threads: get_resources("rseqc_gene_body_coverage")["threads"]
    resources:
        mem_mb=get_resources("rseqc_gene_body_coverage")["mem_mb"],
        runtime=get_resources("rseqc_gene_body_coverage")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/rseqc_gene_body_coverage/{sample}.txt",
    log:
        "results/pipeline_info/logs/rseqc/gene_body_coverage/{sample}.log",
    conda:
        RSEQC_ENV
    shell:
        "geneBody_coverage.py -i {input.aln} -r {input.refgene} "
        "-o {params.prefix} > {log} 2>&1"
