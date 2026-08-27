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


# Only defined when we actually build the index. With build_index=false the
# external index is a plain pre-existing INPUT (STAR_INDEX_DIR points straight
# at it) -- no rule produces that path, so Snakemake never deletes or rewrites
# it, and no multi-GB copy into results/ is needed. See the comment in
# common.smk for why this replaced the old copy-into-results workaround.
if STAR_BUILD_INDEX:

    rule star_index:
        # Builds the STAR genome index once; reused by every sample's
        # alignment.  sjdbOverhang defaults to being auto-detected from the
        # sample sheet's fastq read lengths (see SJDB_OVERHANG in common.smk)
        # unless ref.sjdb_overhang is set to an explicit integer in
        # config.yaml.  Uses the STAR version pinned in
        # config["versions"]["star"] rather than a snakemake-wrapper, so the
        # exact tool version is fully under your control.
        #
        # The fasta/gtf inputs are marked ancient() so an existing index
        # directory is honored as-is instead of being rebuilt whenever the
        # reference files happen to be *newer* than the index.
        input:
            fasta=ancient(FASTA),
            gtf=ancient(GTF),
        output:
            directory(STAR_INDEX_DIR),
        params:
            sjdb_overhang=SJDB_OVERHANG,
            extra=config["star"].get("index_extra", ""),
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
            # If the output directory already contains a valid STAR index
            # (core files SA, SAindex, Genome all present), skip the rebuild.
            # Otherwise build from scratch.  STAR has a long-standing
            # crash-on-exit bug (benign -- output is already written), so a
            # non-zero exit is tolerated if the index files check out.
            "if [ -s {output}/SA ] && [ -s {output}/SAindex ] && "
            "    [ -s {output}/Genome ]; then "
            "  echo 'STAR index already present at {output}; skipping rebuild' "
            "  >> {log}; "
            "else "
            "  mkdir -p {output} && "
            "  (STAR --runMode genomeGenerate "
            "  --genomeDir {output} "
            "  --genomeFastaFiles {input.fasta} "
            "  --sjdbGTFfile {input.gtf} "
            "  --sjdbOverhang {params.sjdb_overhang} "
            "  --runThreadN {threads} "
            "  {params.extra} "
            "  > {log} 2>&1 "
            "  || (echo 'STAR exited non-zero; checking whether the index was "
            "  actually written successfully anyway (see rule comment re: "
            "  known benign STAR exit-time crash)' >> {log}; "
            "  test -s {output}/SA && test -s {output}/SAindex && "
            "  test -s {output}/Genome)); "
            "fi"


rule cleanup_star_index:
    # Removes the STAR genome index after alignment is done, when the user
    # sets outputs.keep_star_index: false.  Saves disk on big runs at the
    # cost of rebuilding the index on the next run.  Runs locally (not
    # submitted to the cluster) -- it is a trivial shell command.
    #
    # Deletes ONLY an index this workflow built itself.  With
    # star.build_index: false the index belongs to the user (it is a plain
    # input, see common.smk), so params.index_dir is empty and the rule is a
    # no-op that just touches its sentinel -- keep_star_index is ignored
    # rather than being allowed to rm -rf someone's shared reference.
    input:
        bams=expand("results/star/{sample}_Aligned.out.bam", sample=SAMPLES),
    output:
        touch("results/pipeline_info/.star_index_cleaned"),
    params:
        keep=KEEP_STAR_INDEX,
        # Empty unless we built it -- the structural guarantee that an
        # external index can never reach the rm -rf below.
        index_dir=STAR_INDEX_DIR if STAR_BUILD_INDEX else "",
    shell:
        "if [ {params.keep} = False ] && [ -n '{params.index_dir}' ] && "
        "    [ -d '{params.index_dir}' ]; then "
        "  echo 'Removing STAR index at {params.index_dir}' && "
        "  rm -rf '{params.index_dir}'; "
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


if CHIMERA_JUNCTION_ENABLED or CHIMERA_ASSEMBLY_ENABLED:

    rule annotation_to_bed:
        # Converts the gene GTF + the curated TE GTF into the BED tracks
        # both chimera screens' breakpoint/exon-overlap tests run against
        # (genes.bed, exons.bed, te.bed) -- shared by chimera_junction.smk
        # and chimera_assembly.smk, so it lives here (built whenever either
        # is enabled) rather than in either one specifically.
        # Pure-python (annotation_to_bed.py), so it runs in the base
        # environment.
        input:
            gtf=GTF,
            te_gtf=TE_GTF,
        output:
            genes="results/reference/genes.bed",
            exons="results/reference/exons.bed",
            te="results/reference/te.bed",
        params:
            outdir="results/reference",
        threads: get_resources("annotation_to_bed")["threads"]
        resources:
            mem_mb=get_resources("annotation_to_bed")["mem_mb"],
            runtime=get_resources("annotation_to_bed")["runtime"],
        benchmark:
            "results/pipeline_info/benchmarks/annotation_to_bed/annotation_to_bed.txt",
        log:
            "results/pipeline_info/logs/reference/annotation_to_bed.log",
        shell:
            "python3 {SCRIPTS_DIR}/annotation_to_bed.py "
            "--gtf {input.gtf} --te-gtf {input.te_gtf} "
            "--outdir {params.outdir} > {log} 2>&1"
