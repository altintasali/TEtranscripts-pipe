rule telocal_locind:
    # Build the TElocal locus index from the TE GTF.  Only triggered when
    # the user does not provide a pre-built .locInd via config
    # telocal.locind -- the telocal rule selects the auto-built path when
    # locind is empty.  Lives under results/telocal/ next to the count
    # tables it feeds; deleted after the last TElocal run when the user sets
    # outputs.keep_telocal_index: false (cleanup_telocal_index below).
    input:
        te_gtf=TE_GTF,
    output:
        # Must end in .locInd -- TElocal rejects any --TE file whose path
        # does not have that suffix.
        "results/telocal/telocal.locInd",
    params:
        # "fast" (default): our reimplementation of TEfeatures.build(),
        # verified against the original -- see build_telocal_index.py's
        # module docstring. "legacy": TElocal_Toolkit's own unmodified
        # build(), an escape hatch (config telocal.indexer).
        indexer=config["telocal"].get("indexer", "fast"),
    threads: get_resources("telocal_locind")["threads"]
    resources:
        mem_mb=get_resources("telocal_locind")["mem_mb"],
        runtime=get_resources("telocal_locind")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_locind/locind.txt",
    log:
        "results/pipeline_info/logs/telocal/locind.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "python3 {SCRIPTS_DIR}/build_telocal_index.py "
        "--gtf {input.te_gtf} "
        "--indexer {params.indexer} "
        "--out {output} > {log} 2>&1"


rule telocal:
    # Per-sample locus-level TE quantification (TElocal). Complements TEcount's
    # subfamily-level quantification by resolving TEs per genomic instance.
    # Uses the same unsorted BAM as TEcount; no re-alignment needed.
    # TElocal has no --outdir flag, so we cd into the output directory before
    # running it; {input.bam} resolves to an absolute path so the cd is safe.
    input:
        bam="results/star/{sample}_Aligned.out.bam",
        gtf=GTF,
        locind=_telocal_locind_path(),
        strandedness=strandedness_input,
    output:
        "results/telocal/{sample}.cntTable.gz",
    params:
        stranded=get_strandedness_param,
        extra=config["telocal"]["extra"],
    threads: get_resources("telocal")["threads"]
    resources:
        mem_mb=get_resources("telocal")["mem_mb"],
        runtime=get_resources("telocal")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal/{sample}.txt",
    log:
        "results/pipeline_info/logs/telocal/{sample}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "( mkdir -p results/telocal &&"
        " TElocal --sortByPos -b {input.bam}"
        " --GTF {input.gtf} --TE {input.locind}"
        " --stranded {params.stranded}"
        " --project {wildcards.sample}"
        " {params.extra}"
        " && mv {wildcards.sample}.cntTable results/telocal/"
        " && gzip -f results/telocal/{wildcards.sample}.cntTable )"
        " > {log} 2>&1"


def telocal_counts_input():
    return [
        f"results/telocal/{s}.cntTable.gz"
        for s in SAMPLES
    ]


rule telocal_summary:
    # Per-sample TElocal summary stats for the MultiQC report (custom
    # content): gene-vs-TE assignment and TE class composition as counts and
    # percentages (telocal_summary_mqc.py). Uses the RAW cntTables.
    input:
        tables=telocal_counts_input(),
    output:
        assignment="results/telocal/qc/telocal_assignment_mqc.json",
        te_class="results/telocal/qc/telocal_te_class_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("telocal")["threads"]
    resources:
        mem_mb=get_resources("telocal")["mem_mb"],
        runtime=get_resources("telocal")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_summary/telocal_summary.txt",
    log:
        "results/pipeline_info/logs/telocal/summary.log",
    shell:
        "python3 {SCRIPTS_DIR}/telocal_summary_mqc.py "
        "--tables {input.tables} "
        "--samples {params.samples} "
        "--out-assignment {output.assignment} "
        "--out-class {output.te_class} > {log} 2>&1"


rule cleanup_telocal_index:
    # Removes the auto-built TElocal .locInd index after every sample's
    # TElocal run finished, when the user sets outputs.keep_telocal_index:
    # false.  Saves disk at the cost of re-building the index on the next
    # run.  A user-provided telocal.locind path is never touched -- it is
    # only ever read.  Runs locally (not submitted to the cluster) -- it is
    # a trivial shell command.
    input:
        tables=telocal_counts_input(),
    output:
        touch("results/pipeline_info/.telocal_index_cleaned"),
    params:
        keep=KEEP_TELOCAL_INDEX,
        locind=_telocal_locind_path(),
        auto_built=not _telocal_locind_cfg,
    shell:
        "if [ {params.keep} = False ] && [ {params.auto_built} = True ] && "
        "[ -f {params.locind} ]; then "
        "  echo 'Removing auto-built TElocal index at {params.locind}' && "
        "  rm -f {params.locind}; "
        "fi"


rule telocal_counts:
    # Merges every sample's TElocal count table into the locus x sample
    # counts matrix (tecount_counts.py --key-style telocal) that feeds the
    # sample-QC stage.  --feature-class restricts the matrix to TE loci
    # (default), genes, or all features; only this QC-view matrix is
    # filtered, never the per-sample cntTables.
    input:
        tables=telocal_counts_input(),
    output:
        counts="results/telocal/counts_matrix.tsv.gz",
    params:
        sample_names=lambda wc, input: " ".join(SAMPLES),
        feature_class=TELOCAL_QC["feature_class"],
    threads: get_resources("telocal_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("telocal_counts"),
        runtime=get_resources("telocal_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_counts/telocal_counts.txt",
    log:
        "results/pipeline_info/logs/telocal/counts.log",
    shell:
        "python3 {SCRIPTS_DIR}/tecount_counts.py "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--key-style telocal "
        "--feature-class {params.feature_class} "
        "--out-counts {output.counts} > {log} 2>&1"


rule telocal_qc_transform:
    # Normalizes the TElocal counts matrix for the QC view (sample_qc.R
    # --transform telocal).  Defaults to log2 (locus-level matrices are far
    # too large for vst/rlog); vst/rlog remain selectable.  Filters in
    # telocal.qc apply only to this view.
    input:
        counts="results/telocal/counts_matrix.tsv.gz",
    output:
        "results/telocal/qc/{transform}_counts.tsv.gz",
    params:
        samples=config["samples"],
        min_samples_present=TELOCAL_QC["min_samples_present"],
        min_total_counts=TELOCAL_QC["min_total_counts"],
    threads: get_resources("telocal_qc_transform")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("telocal_qc_transform"),
        runtime=get_resources("telocal_qc_transform")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_qc_transform/{transform}.txt",
    log:
        "results/pipeline_info/logs/telocal/qc/transform_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/sample_qc.R "
        "--transform telocal {input.counts} {params.samples} {wildcards.transform} "
        "{params.min_samples_present} {params.min_total_counts} "
        "{output} > {log} 2>&1"


rule telocal_qc:
    # PCA scatter + sample-to-sample distance heatmap of the transformed
    # TElocal counts, colored by condition (sample sheet's "condition"
    # column; absent -> one "all" group), emitted as MultiQC custom-content
    # JSON (ids telocal_sample_qc_pca / telocal_sample_qc_heatmap, ordered
    # inside the custom_content module by multiqc_config.yaml).
    input:
        transformed="results/telocal/qc/{transform}_counts.tsv.gz",
    output:
        pca="results/telocal/qc/pca_{transform}_mqc.json",
        heatmap="results/telocal/qc/heatmap_{transform}_mqc.json",
    params:
        samples=config["samples"],
        min_events=TELOCAL_QC["min_events"],
    threads: get_resources("telocal_qc")["threads"]
    resources:
        mem_mb=get_resources("telocal_qc")["mem_mb"],
        runtime=get_resources("telocal_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/telocal_qc/{transform}.txt",
    log:
        "results/pipeline_info/logs/telocal/qc/plots_{transform}.log",
    conda:
        TETRANSCRIPTS_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/sample_qc.R "
        "--plots telocal {input.transformed} {params.samples} {params.min_events} "
        "{wildcards.transform} {output.pca} {output.heatmap} > {log} 2>&1"
