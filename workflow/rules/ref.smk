rule gunzip_reference:
    # Decompresses a gzipped reference file (fasta/gtf/te_gtf) so downstream
    # tools that expect plain text (STAR's --genomeFastaFiles/--sjdbGTFfile,
    # TEcount/TEtranscripts' --GTF/--TE) can consume it directly. Only
    # triggered for whichever of fasta/gtf/te_gtf are actually given as
    # .gz in config.yaml -- see REFERENCE_GZ_SOURCES in common.smk. Outputs
    # go to ref.decompressed_dir (default: results/pipeline_info/
    # ref_decompressed -- shared storage, unlike a node-local /tmp which
    # breaks cluster runs) and are temp()-wrapped so they are deleted once
    # star_index/gtf_to_genepred/tecount are done with them.
    input:
        lambda wc: REFERENCE_GZ_SOURCES[wc.stem],
    output:
        temp(f"{DECOMPRESS_DIR}/{{stem}}"),
    threads: get_resources("gunzip_reference")["threads"]
    resources:
        mem_mb=get_resources("gunzip_reference")["mem_mb"],
        runtime=get_resources("gunzip_reference")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/gunzip_reference/{stem}.txt",
    log:
        "results/pipeline_info/logs/gunzip/{stem}.log",
    shell:
        "mkdir -p {DECOMPRESS_DIR} && "
        "gunzip -c {input} > {output} 2> {log}"


rule star_index:
    # Builds the STAR genome index once; reused by every sample's alignment.
    # sjdbOverhang defaults to being auto-detected from the sample sheet's
    # fastq read lengths (see SJDB_OVERHANG in common.smk) unless
    # ref.sjdb_overhang is set to an explicit integer in config.yaml.
    # Uses the STAR version pinned in config["versions"]["star"] (see the
    # generated env in common.smk) rather than a snakemake-wrapper, so the
    # exact tool version is fully under your control.
    #
    # The fasta/gtf inputs are marked ancient() when build_index=true: an
    # existing index directory is honored as-is instead of being rebuilt
    # whenever the reference files happen to be *newer* than the index
    # (Snakemake's up-to-date check compares mtimes, which would otherwise
    # re-trigger the index build every run for a shared/prebuilt index --
    # see the "STAR index" README note).  When build_index=false the inputs
    # are empty lists so no decompression (gunzip_reference) is triggered
    # and the output directory already exists, so Snakemake skips the rule
    # entirely.
    input:
        fasta=ancient(FASTA) if STAR_BUILD_INDEX else [],
        gtf=ancient(GTF) if STAR_BUILD_INDEX else [],
    output:
        directory(STAR_INDEX_DIR),
    params:
        sjdb_overhang=SJDB_OVERHANG,
        extra=config["star"].get("index_extra", ""),
        build_index=STAR_BUILD_INDEX,
    threads: get_resources("star_index")["threads"]
    resources:
        mem_mb=get_resources("star_index")["mem_mb"],
        runtime=get_resources("star_index")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/star_index/star_index.txt",
    log:
        "results/pipeline_info/logs/star/index.log",
    conda:
        STAR_ENV
    shell:
        # When build_index=true: if the output directory already contains a
        # valid STAR index (core files SA, SAindex, Genome all present), skip
        # the rebuild entirely. Otherwise, build from scratch. STAR has a
        # long-standing, still-unresolved crash-on-exit bug ("double free or
        # corruption" / segfault while freeing memory *after* all output has
        # already been written and closed) -- confirmed benign by STAR's own
        # author: https://groups.google.com/g/rna-star/c/3_ckDieghws ("this
        # problem is happening after STAR finished all calculations, so it
        # does not affect the results"). So if STAR exits non-zero, don't
        # trust that alone -- check whether the core index files actually
        # landed on disk before treating it as a real failure.
        #
        # When build_index=false: the output directory already exists and
        # the rule is skipped by Snakemake before the shell runs, but if
        # --rerun-incomplete forces re-execution this branch is a safe no-op.
        "if [ {params.build_index} = True ]; then "
        "if [ -s {output}/SA ] && [ -s {output}/SAindex ] && "
        "[ -s {output}/Genome ]; then "
        "echo 'STAR index already present at {output}; skipping rebuild' "
        ">> {log}; "
        "else "
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
        "test -s {output}/Genome)); "
        "fi; "
        "else "
        "echo 'star.build_index=false; index already exists' >> {log} || true; "
        "fi"


rule gtf_to_genepred:
    # RSeQC's infer_experiment.py needs a BED12 reference model, not a GTF.
    # UCSC's gtfToGenePred + genePredToBed is the standard conversion route.
    input:
        gtf=GTF,
    output:
        genepred=temp("results/reference/annotation.genePred"),
    threads: get_resources("gtf_to_genepred")["threads"]
    resources:
        mem_mb=get_resources("gtf_to_genepred")["mem_mb"],
        runtime=get_resources("gtf_to_genepred")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/gtf_to_genepred/gtf_to_genepred.txt",
    log:
        "results/pipeline_info/logs/rseqc/gtf_to_genepred.log",
    conda:
        UCSC_TOOLS_ENV
    shell:
        "gtfToGenePred -genePredExt -ignoreGroupsWithoutExons "
        "{input.gtf} {output.genepred} > {log} 2>&1"


rule genepred_to_bed12:
    input:
        genepred="results/reference/annotation.genePred",
    output:
        bed12="results/reference/annotation.bed12",
    threads: get_resources("genepred_to_bed12")["threads"]
    resources:
        mem_mb=get_resources("genepred_to_bed12")["mem_mb"],
        runtime=get_resources("genepred_to_bed12")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/genepred_to_bed12/genepred_to_bed12.txt",
    log:
        "results/pipeline_info/logs/rseqc/genepred_to_bed12.log",
    conda:
        UCSC_TOOLS_ENV
    shell:
        "genePredToBed {input.genepred} {output.bed12} > {log} 2>&1"
