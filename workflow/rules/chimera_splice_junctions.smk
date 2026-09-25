# -----------------------------------------------------------------------------
# Chimera SJ.out.tab junction screen: gene-TE chimera detection from STAR's
# own normal splice junctions -- a third, independent evidence source
# alongside chimera_chimeric_reads.smk (STAR chimeric-junction reads) and
# chimera_assembly.smk (StringTie assembly structure).
#
# OFF BY DEFAULT: brand new, unvalidated on real data. STAR only writes a
# chimeric-junction record when a read can't be explained by one linear
# (possibly spliced) alignment -- a TE that splices into a gene via an
# ordinary, canonical, nearby intron aligns as a completely normal spliced
# read and never reaches chimera_chimeric_reads at all. chimera_assembly.smk exists
# to catch that same blind spot from StringTie's assembled transcript
# structure, which needs enough per-sample coverage to successfully
# assemble a multi-exon transcript; this screen catches it instead at the
# individual read-junction level, straight from results/star/{sample}_
# SJ.out.tab -- already an UNCONDITIONAL output of the main star_align rule
# (align.smk) for every sample, so no extra STAR pass and no assembly step
# are needed here.
#
# genes.bed/exons.bed/te.bed are built by ref.smk's annotation_to_bed rule
# (shared with chimera_chimeric_reads.smk/chimera_assembly.smk).
#
# Rules:
#   chimera_splice_junctions_classify      per-sample SJ.out.tab annotation (+ the gene-TE subset)
#   chimera_splice_junctions_counts        merge per-sample tables -> all_events + counts
#                            + te-gene-junctions
#   chimera_splice_junctions_qc_transform  count-matrix transform for the sample-QC view
#   chimera_splice_junctions_qc            sample-QC PCA / sample-distance plots
# -----------------------------------------------------------------------------
import os

WRITE_SPLICE_JUNCTIONS_COUNTS = bool(config["chimera"]["splice_junctions"]["outputs"]["write_counts_matrix"])
CHIMERA_SPLICE_JUNCTIONS_QC = config["chimera"]["splice_junctions"]["qc"]


def chimera_splice_junctions_counts_input():
    return [
        f"results/chimera/splice_junctions/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]


def all_chimera_splice_junctions_outputs():
    """Chimera-SJ-junction-screen artifacts for the `all` target (Snakefile)."""
    files = [
        f"results/chimera/splice_junctions/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        f"results/chimera/splice_junctions/per_sample/{s}_junctions_te-gene-junctions.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/splice_junctions/all_events.tsv.gz",
        "results/chimera/splice_junctions/te-gene-junctions.tsv.gz",
    ]
    if WRITE_SPLICE_JUNCTIONS_COUNTS:
        files += [
            "results/chimera/splice_junctions/counts_matrix.tsv.gz",
            "results/chimera/splice_junctions/cpm_matrix.tsv.gz",
        ]
    return files


def all_chimera_splice_junctions_sample_qc_outputs():
    """Sample-QC artifacts for the `all` target (Snakefile). Only produced
    when a counts matrix is written (it is the QC view's input).

    Written under the SHARED results/chimera/qc/ directory, splice_junctions_-
    prefixed, matching chimera_assembly.smk's convention -- MultiQC's search
    directories (multiqc_config.yaml) only scan results/chimera/qc, not a
    per-screen qc/ subdirectory, so anything written elsewhere silently never
    reaches the report."""
    if not WRITE_SPLICE_JUNCTIONS_COUNTS:
        return []
    transform = CHIMERA_SPLICE_JUNCTIONS_QC["pca_transform"]
    return [
        f"results/chimera/qc/splice_junctions_{transform}_counts.tsv.gz",
        f"results/chimera/qc/splice_junctions_pca_{transform}_mqc.json",
        f"results/chimera/qc/splice_junctions_heatmap_{transform}_mqc.json",
    ]


rule chimera_splice_junctions_classify:
    # Annotates one sample's STAR normal splice junctions (SJ.out.tab)
    # against the gene/TE tracks. See classify_chimera_splice_junctions.py for the full
    # column spec. Pure-python, so it runs in the base environment.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/classify_chimera_splice_junctions.py",
        sj="results/star/{sample}_SJ.out.tab",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
        strandedness=strandedness_input,
    output:
        junctions="results/chimera/splice_junctions/per_sample/{sample}_junctions.tsv.gz",
        te_gene_junctions="results/chimera/splice_junctions/per_sample/{sample}_junctions_te-gene-junctions.tsv.gz",
    params:
        tolerance=config["chimera"]["splice_junctions"]["breakpoint_tolerance"],
        min_unique_reads=config["chimera"]["splice_junctions"]["min_unique_reads"],
        canonical_flag=(
            "--require-canonical"
            if config["chimera"]["splice_junctions"]["require_canonical"] else ""
        ),
        library=get_strandedness_param,
    threads: get_resources("chimera_splice_junctions_classify")["threads"]
    resources:
        mem_mb=get_resources("chimera_splice_junctions_classify")["mem_mb"],
        runtime=get_resources("chimera_splice_junctions_classify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_splice_junctions_classify/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_splice_junctions/classify/{sample}.log",
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


rule chimera_splice_junctions_counts:
    # Merges every sample's junction table into the all-events catalog, the
    # event x sample raw-counts matrix (STAR's own unique_reads per
    # junction), and its CPM-normalized sibling (chimera_splice_junctions_counts.py).
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_splice_junctions_counts.py",
        tables=chimera_splice_junctions_counts_input(),
    output:
        events="results/chimera/splice_junctions/all_events.tsv.gz",
        te_events="results/chimera/splice_junctions/te-gene-junctions.tsv.gz",
        **({"counts": "results/chimera/splice_junctions/counts_matrix.tsv.gz",
            "cpm": "results/chimera/splice_junctions/cpm_matrix.tsv.gz"}
           if WRITE_SPLICE_JUNCTIONS_COUNTS else {}),
    params:
        sample_names=" ".join(SAMPLES),
        counts_flag=(
            "--out-counts results/chimera/splice_junctions/counts_matrix.tsv.gz "
            "--out-cpm results/chimera/splice_junctions/cpm_matrix.tsv.gz"
            if WRITE_SPLICE_JUNCTIONS_COUNTS else ""
        ),
    threads: get_resources("chimera_splice_junctions_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_splice_junctions_counts"),
        runtime=get_resources("chimera_splice_junctions_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_splice_junctions_counts/chimera_splice_junctions_counts.txt",
    log:
        "results/pipeline_info/logs/chimera_splice_junctions/counts.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--out-events {output.events} "
        "--out-te-events {output.te_events} "
        "{params.counts_flag} > {log} 2>&1"


if WRITE_SPLICE_JUNCTIONS_COUNTS:

    rule chimera_splice_junctions_qc_transform:
        input:
            # Declared so that EDITING the script re-runs the rule.
            # Snakemake's code trigger hashes the shell command STRING,
            # not the file it names, so without this an edit to the
            # script leaves stale outputs in place silently.
            script=f"{SCRIPTS_DIR}/sample_qc.R",
            counts="results/chimera/splice_junctions/counts_matrix.tsv.gz",
        output:
            "results/chimera/qc/splice_junctions_{transform}_counts.tsv.gz",
        params:
            samples=config["samples"],
            min_samples_present=CHIMERA_SPLICE_JUNCTIONS_QC["min_samples_present"],
            min_total_counts=CHIMERA_SPLICE_JUNCTIONS_QC["min_total_counts"],
        threads: get_resources("chimera_splice_junctions_qc_transform")["threads"]
        resources:
            mem_mb=get_scaled_mem_mb("chimera_splice_junctions_qc_transform"),
            runtime=get_resources("chimera_splice_junctions_qc_transform")["runtime"],
        benchmark:
            "results/pipeline_info/benchmarks/chimera_splice_junctions_qc_transform/{transform}.txt",
        log:
            "results/pipeline_info/logs/chimera_splice_junctions/qc_transform_{transform}.log",
        conda:
            TETRANSCRIPTS_ENV
        shell:
            "Rscript {input.script} "
            "--transform splice_junctions {input.counts} {params.samples} {wildcards.transform} "
            "{params.min_samples_present} {params.min_total_counts} "
            "{output} > {log} 2>&1"


    rule chimera_splice_junctions_qc:
        input:
            # Declared so that EDITING the script re-runs the rule.
            # Snakemake's code trigger hashes the shell command STRING,
            # not the file it names, so without this an edit to the
            # script leaves stale outputs in place silently.
            script=f"{SCRIPTS_DIR}/sample_qc.R",
            transformed="results/chimera/qc/splice_junctions_{transform}_counts.tsv.gz",
        output:
            pca="results/chimera/qc/splice_junctions_pca_{transform}_mqc.json",
            heatmap="results/chimera/qc/splice_junctions_heatmap_{transform}_mqc.json",
        params:
            samples=config["samples"],
            min_events=CHIMERA_SPLICE_JUNCTIONS_QC["min_events"],
        threads: get_resources("chimera_splice_junctions_qc")["threads"]
        resources:
            mem_mb=get_scaled_mem_mb("chimera_splice_junctions_qc"),
            runtime=get_resources("chimera_splice_junctions_qc")["runtime"],
        benchmark:
            "results/pipeline_info/benchmarks/chimera_splice_junctions_qc/{transform}.txt",
        log:
            "results/pipeline_info/logs/chimera_splice_junctions/qc_plots_{transform}.log",
        conda:
            TETRANSCRIPTS_ENV
        shell:
            "Rscript {input.script} "
            "--plots splice_junctions {input.transformed} {params.samples} {params.min_events} "
            "{wildcards.transform} {output.pca} {output.heatmap} > {log} 2>&1"
