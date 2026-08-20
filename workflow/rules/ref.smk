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
    # build_index=true:  the fasta/gtf inputs are marked ancient() so an
    #   existing index directory is honored as-is instead of being rebuilt
    #   whenever the reference files happen to be *newer* than the index.
    #
    # build_index=false: the user provides a pre-built external index.  The
    #   inputs are empty (no gunzip_reference triggered) and the external
    #   index is copied into the output directory (idempotent -- skipped if
    #   already present).  The copy is persistent (not temp()-wrapped) to
    #   avoid accidental deletion of the original external index -- see the
    #   validation in common.smk that prevents source/destination overlap.
    input:
        fasta=ancient(FASTA) if STAR_BUILD_INDEX else [],
        gtf=ancient(GTF) if STAR_BUILD_INDEX else [],
    output:
        directory(STAR_INDEX_DIR),
    params:
        sjdb_overhang=SJDB_OVERHANG,
        extra=config["star"].get("index_extra", ""),
        build_index=STAR_BUILD_INDEX,
        source_index=STAR_INDEX_SOURCE,
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
        # build_index=true:  if the output directory already contains a
        # valid STAR index (core files SA, SAindex, Genome all present), skip
        # the rebuild.  Otherwise build from scratch.  STAR has a long-
        # standing crash-on-exit bug (benign -- output is already written).
        # We tolerate a non-zero exit by checking the index files.
        #
        # build_index=false: copy the external index into results/ (one-time
        # cost).  Idempotent -- if the copy already exists, skip.  The copy
        # is persistent (not temp()-wrapped) so the original index is never
        # at risk of accidental deletion.
        "if [ {params.build_index} = True ]; then "
        "  if [ -s {output}/SA ] && [ -s {output}/SAindex ] && "
        "      [ -s {output}/Genome ]; then "
        "    echo 'STAR index already present at {output}; skipping rebuild' "
        "    >> {log}; "
        "  else "
        "    mkdir -p {output} && "
        "    (STAR --runMode genomeGenerate "
        "    --genomeDir {output} "
        "    --genomeFastaFiles {input.fasta} "
        "    --sjdbGTFfile {input.gtf} "
        "    --sjdbOverhang {params.sjdb_overhang} "
        "    --runThreadN {threads} "
        "    {params.extra} "
        "    > {log} 2>&1 "
        "    || (echo 'STAR exited non-zero; checking whether the index was "
        "    actually written successfully anyway (see rule comment re: known "
        "    benign STAR exit-time crash)' >> {log}; "
        "    test -s {output}/SA && test -s {output}/SAindex && "
        "    test -s {output}/Genome)); "
        "  fi; "
        "else "
        "  mkdir -p {output} && "
        "  if [ -s {output}/SA ] && [ -s {output}/SAindex ] && "
        "      [ -s {output}/Genome ]; then "
        "    echo 'STAR index already present at {output}; skipping copy' "
        "    >> {log}; "
        "  else "
        "    echo 'Copying external STAR index from {params.source_index} "
        "    to {output}' >> {log} && "
        "    cp -a {params.source_index}/* {output}/; "
        "  fi; "
        "fi"


rule cleanup_star_index:
    # Removes the STAR genome index after alignment is done, when the user
    # sets outputs.keep_star_index: false.  Saves disk on big runs at the
    # cost of re-copying (or rebuilding) the index on the next run.
    input:
        bams=expand("results/star/{sample}_Aligned.out.bam", sample=SAMPLES),
    output:
        touch("results/pipeline_info/.star_index_cleaned"),
    params:
        keep=KEEP_STAR_INDEX,
        index_dir=STAR_INDEX_DIR,
    resources:
        mem_mb=100,
        runtime=1,
    shell:
        "if [ {params.keep} = False ] && [ -d {params.index_dir} ]; then "
        "  echo 'Removing STAR index at {params.index_dir}' && "
        "  rm -rf {params.index_dir}; "
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
