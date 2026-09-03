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


# How many transcripts geneBody_coverage.py is given (see the rule below for
# why this matters so much). Not a config key: it sets the precision of a QC
# curve, not an analysis result -- 1000 genes is where the aggregate stops
# being noisy, and RSeQC's own housekeeping-gene BEDs, its documented answer
# to this same slowness, are only ~3.8k transcripts.
GENE_BODY_COVERAGE_TRANSCRIPTS = 1000


rule rseqc_gene_body_coverage:
    # 5'->3' gene body coverage -- flags RNA degradation (3' bias) or
    # incomplete reverse-transcription (5' bias); MultiQC's RSeQC module
    # picks up geneBody_coverage.py's .txt table as its own subsection.
    # geneBody_coverage.py names outputs from -o's prefix (.txt/.r/.pdf);
    # only the .txt table is a tracked output -- the .r/.pdf are RSeQC's own
    # plotting scaffolding, unused since MultiQC parses the .txt directly.
    #
    # The BED12 is thinned to ~GENE_BODY_COVERAGE_TRANSCRIPTS lines first,
    # because handing this tool a full transcriptome model (~250k transcripts
    # for GENCODE) costs hours per sample. From RSeQC 5.0.4's source: per
    # transcript it samples 100 percentile points from the exonic bases, then
    # pileups from the first to the LAST exonic base -- the whole genomic
    # span, introns included -- in a pure-Python loop. So the cost is
    # sum(genomic span x depth) over transcripts, single-threaded (no
    # internal parallelism, so `threads:` cannot help). Worse, its per-gene
    # key includes the coordinates, so every isoform is a separate entry that
    # re-pileups the same locus from scratch.
    #
    # Taking every k-th line handles both: the BED12 is coordinate-sorted, so
    # a stride spreads the picks across the whole genome (unlike `head`,
    # which would take one clustered stretch of chromosome 1) and, since a
    # gene's isoforms are adjacent lines, it lands on at most one isoform per
    # gene. k is derived per run, so the target holds for any annotation
    # size; it floors at 1, so a model already smaller than the target keeps
    # every line and the count only ever errs high (at most ~2x) rather than
    # low. The other two RSeQC rules stream the BAM once and are not slow, so
    # they keep using the full annotation.bed12.
    input:
        aln="results/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        refgene="results/reference/annotation.bed12",
    output:
        "results/rseqc/{sample}.geneBodyCoverage.txt",
    params:
        prefix="results/rseqc/{sample}",
        n_transcripts=GENE_BODY_COVERAGE_TRANSCRIPTS,
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
        "bed={resources.tmpdir}/{wildcards.sample}.genebody.bed12; "
        "k=$(( $(wc -l < {input.refgene}) / {params.n_transcripts} )); "
        "k=$(( k > 0 ? k : 1 )); "
        "awk -v k=\"$k\" 'NR % k == 0' {input.refgene} > \"$bed\" && "
        "geneBody_coverage.py -i {input.aln} -r \"$bed\" "
        "-o {params.prefix} > {log} 2>&1"
