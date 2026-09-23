# -----------------------------------------------------------------------------
# Chimera SJ.out.tab junction screen: gene-TE chimera detection from STAR's
# own normal splice junctions -- a third, independent evidence source
# alongside chimera_reads.smk (STAR chimeric-junction reads) and
# chimera_assembly.smk (StringTie assembly structure).
#
# OFF BY DEFAULT: brand new, unvalidated on real data. STAR only writes a
# chimeric-junction record when a read can't be explained by one linear
# (possibly spliced) alignment -- a TE that splices into a gene via an
# ordinary, canonical, nearby intron aligns as a completely normal spliced
# read and never reaches chimera_reads at all. chimera_assembly.smk exists
# to catch that same blind spot from StringTie's assembled transcript
# structure, which needs enough per-sample coverage to successfully
# assemble a multi-exon transcript; this screen catches it instead at the
# individual read-junction level, straight from results/star/{sample}_
# SJ.out.tab -- already an UNCONDITIONAL output of the main star_align rule
# (align.smk) for every sample, so no extra STAR pass and no assembly step
# are needed here.
#
# genes.bed/exons.bed/te.bed are built by ref.smk's annotation_to_bed rule
# (shared with chimera_reads.smk/chimera_assembly.smk).
#
# Rules:
#   chimera_sj_classify      per-sample SJ.out.tab annotation (+ the gene-TE subset)
#   chimera_sj_counts        merge per-sample tables -> all_events + counts
#                            + te-gene-junctions
#   chimera_sj_qc_transform  count-matrix transform for the sample-QC view
#   chimera_sj_qc            sample-QC PCA / sample-distance plots
# -----------------------------------------------------------------------------
import os

WRITE_SJ_COUNTS = bool(config["chimera"]["sj_junctions"]["outputs"]["write_counts_matrix"])
CHIMERA_SJ_QC = config["chimera"]["sj_junctions"]["qc"]


def chimera_sj_counts_input():
    return [
        f"results/chimera/sj/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]


def all_chimera_sj_outputs():
    """Chimera-SJ-junction-screen artifacts for the `all` target (Snakefile)."""
    files = [
        f"results/chimera/sj/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        f"results/chimera/sj/per_sample/{s}_junctions_te-gene-junctions.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/sj/all_events.tsv.gz",
        "results/chimera/sj/te-gene-junctions.tsv.gz",
    ]
    if WRITE_SJ_COUNTS:
        files += [
            "results/chimera/sj/counts_matrix.tsv.gz",
            "results/chimera/sj/cpm_matrix.tsv.gz",
        ]
    return files


def all_chimera_sj_sample_qc_outputs():
    """Sample-QC artifacts for the `all` target (Snakefile). Only produced
    when a counts matrix is written (it is the QC view's input)."""
    if not WRITE_SJ_COUNTS:
        return []
    transform = CHIMERA_SJ_QC["pca_transform"]
    return [
        f"results/chimera/sj/qc/{transform}_counts.tsv.gz",
        f"results/chimera/sj/qc/pca_{transform}_mqc.json",
        f"results/chimera/sj/qc/heatmap_{transform}_mqc.json",
    ]


rule chimera_sj_classify:
    # Annotates one sample's STAR normal splice junctions (SJ.out.tab)
    # against the gene/TE tracks. See classify_chimera_sj.py for the full
    # column spec. Pure-python, so it runs in the base environment.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/classify_chimera_sj.py",
        sj="results/star/{sample}_SJ.out.tab",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
        strandedness=strandedness_input,
    output:
        junctions="results/chimera/sj/per_sample/{sample}_junctions.tsv.gz",
        te_gene_junctions="results/chimera/sj/per_sample/{sample}_junctions_te-gene-junctions.tsv.gz",
    params:
        tolerance=config["chimera"]["sj_junctions"]["breakpoint_tolerance"],
        min_unique_reads=config["chimera"]["sj_junctions"]["min_unique_reads"],
        canonical_flag=(
            "--require-canonical"
            if config["chimera"]["sj_junctions"]["require_canonical"] else ""
        ),
        library=get_strandedness_param,
    threads: get_resources("chimera_sj_classify")["threads"]
    resources:
        mem_mb=get_resources("chimera_sj_classify")["mem_mb"],
        runtime=get_resources("chimera_sj_classify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_sj_classify/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_sj/classify/{sample}.log",
    shell:
        "python3 {input.script} "
        "--sj {input.sj} "
        "--genes {input.genes} --exons {input.exons} --te {input.te} "
        "--sample {wildcards.sample} "
        "--breakpoint-tolerance {params.tolerance} "
        "--min-unique-reads {params.min_unique_reads} "
        "{params.canonical_flag} "
        "--library-strandedness {params.library} "
        "--out {output.junctions} "
        "--te-out {output.te_gene_junctions} > {log} 2>&1"


rule chimera_sj_counts:
    # Merges every sample's junction table into the all-events catalog, the
    # event x sample raw-counts matrix (STAR's own unique_reads per
    # junction), and its CPM-normalized sibling (chimera_sj_counts.py).
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_sj_counts.py",
        tables=chimera_sj_counts_input(),
    output:
        events="results/chimera/sj/all_events.tsv.gz",
        te_events="results/chimera/sj/te-gene-junctions.tsv.gz",
        **({"counts": "results/chimera/sj/counts_matrix.tsv.gz",
            "cpm": "results/chimera/sj/cpm_matrix.tsv.gz"}
           if WRITE_SJ_COUNTS else {}),
    params:
        sample_names=" ".join(SAMPLES),
        counts_flag=(
            "--out-counts results/chimera/sj/counts_matrix.tsv.gz "
            "--out-cpm results/chimera/sj/cpm_matrix.tsv.gz"
            if WRITE_SJ_COUNTS else ""
        ),
    threads: get_resources("chimera_sj_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_sj_counts"),
        runtime=get_resources("chimera_sj_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_sj_counts/chimera_sj_counts.txt",
    log:
        "results/pipeline_info/logs/chimera_sj/counts.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--out-events {output.events} "
        "--out-te-events {output.te_events} "
        "{params.counts_flag} > {log} 2>&1"


if WRITE_SJ_COUNTS:

    rule chimera_sj_qc_transform:
        input:
            # Declared so that EDITING the script re-runs the rule.
            # Snakemake's code trigger hashes the shell command STRING,
            # not the file it names, so without this an edit to the
            # script leaves stale outputs in place silently.
            script=f"{SCRIPTS_DIR}/sample_qc.R",
            counts="results/chimera/sj/counts_matrix.tsv.gz",
        output:
            "results/chimera/sj/qc/{transform}_counts.tsv.gz",
        params:
            samples=config["samples"],
            min_samples_present=CHIMERA_SJ_QC["min_samples_present"],
            min_total_counts=CHIMERA_SJ_QC["min_total_counts"],
        threads: get_resources("chimera_sj_qc_transform")["threads"]
        resources:
            mem_mb=get_scaled_mem_mb("chimera_sj_qc_transform"),
            runtime=get_resources("chimera_sj_qc_transform")["runtime"],
        benchmark:
            "results/pipeline_info/benchmarks/chimera_sj_qc_transform/{transform}.txt",
        log:
            "results/pipeline_info/logs/chimera_sj/qc_transform_{transform}.log",
        conda:
            TETRANSCRIPTS_ENV
        shell:
            "Rscript {input.script} "
            "--transform sj {input.counts} {params.samples} {wildcards.transform} "
            "{params.min_samples_present} {params.min_total_counts} "
            "{output} > {log} 2>&1"


    rule chimera_sj_qc:
        input:
            # Declared so that EDITING the script re-runs the rule.
            # Snakemake's code trigger hashes the shell command STRING,
            # not the file it names, so without this an edit to the
            # script leaves stale outputs in place silently.
            script=f"{SCRIPTS_DIR}/sample_qc.R",
            transformed="results/chimera/sj/qc/{transform}_counts.tsv.gz",
        output:
            pca="results/chimera/sj/qc/pca_{transform}_mqc.json",
            heatmap="results/chimera/sj/qc/heatmap_{transform}_mqc.json",
        params:
            samples=config["samples"],
            min_events=CHIMERA_SJ_QC["min_events"],
        threads: get_resources("chimera_sj_qc")["threads"]
        resources:
            mem_mb=get_scaled_mem_mb("chimera_sj_qc"),
            runtime=get_resources("chimera_sj_qc")["runtime"],
        benchmark:
            "results/pipeline_info/benchmarks/chimera_sj_qc/{transform}.txt",
        log:
            "results/pipeline_info/logs/chimera_sj/qc_plots_{transform}.log",
        conda:
            TETRANSCRIPTS_ENV
        shell:
            "Rscript {input.script} "
            "--plots sj {input.transformed} {params.samples} {params.min_events} "
            "{wildcards.transform} {output.pca} {output.heatmap} > {log} 2>&1"
