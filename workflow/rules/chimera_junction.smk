# -----------------------------------------------------------------------------
# Chimera screen: gene-TE chimeric junction detection and annotation.
#
# EXPERIMENTAL: newer addition to the pipeline -- junction classification and
# the sample-QC view may still change between releases. Validate output before
# relying on it for published results.
#
# The STAR alignment run already emits {sample}_Chimeric.out.junction (see
# align.smk --chimOutType). This stage annotates those junctions against the
# gene/TE tracks and aggregates them into the all-events catalog and the
# counts matrix the sample-QC stage (sample_qc.smk) consumes.
#
# The genes.bed/exons.bed/te.bed tracks this stage annotates against are
# built by ref.smk's annotation_to_bed rule (shared with chimera_assembly.smk,
# so it lives in ref.smk rather than here -- it must exist whenever EITHER
# chimera screen is enabled, not just this one).
#
# Per-sample tables are two row-sets (all junctions / the gene-TE chimera
# subset) x two column-sets (without / with TElocal columns):
#   {sample}_junctions.tsv.gz                              base
#   {sample}_junctions_te-gene-chimeras.tsv.gz             fewer rows
#   {sample}_junctions_with-telocal.tsv.gz                 more columns
#   {sample}_junctions_te-gene-chimeras_with-telocal.tsv.gz  both
# (the two _with-telocal ones only when telocal is enabled).
#
# Rules:
#   parse_chimeric_junctions per-sample junction annotation (+ the gene-TE subset)
#   chimera_counts           merge per-sample tables -> all_events + counts
#                            + te-gene-chimeras
#   junction_qc              per-sample QC summary (MultiQC custom content)
#   chimera_igv_bed          per-sample IGV track (config-gated)
# -----------------------------------------------------------------------------
import os

WRITE_COUNTS = bool(config["chimera"]["junction"]["outputs"]["write_counts_matrix"])
WRITE_IGV_BED = bool(config["chimera"]["junction"]["outputs"]["write_igv_bed"])
# CHIMERA_QC now lives in common/runtime.smk -- the assembly screen needs it
# too, and this file is not included when the junction screen is off.


def chimera_counts_input():
    if TELOCAL_ENABLED:
        return [
            f"results/chimera/read_evidence/per_sample/{s}_junctions_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
    return [
        f"results/chimera/read_evidence/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]


def all_chimera_outputs():
    """Chimera artifacts for the `all` target (Snakefile)."""
    files = [
        f"results/chimera/read_evidence/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]
    if TELOCAL_ENABLED:
        files += [
            f"results/chimera/read_evidence/per_sample/{s}_junctions_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
        files += [
            f"results/chimera/read_evidence/per_sample/{s}_junctions_te-gene-chimeras_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
    files += [
        f"results/chimera/read_evidence/per_sample/{s}_junctions_te-gene-chimeras.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/read_evidence/all_events.tsv.gz",
        "results/chimera/read_evidence/te-gene-chimeras.tsv.gz",
        "results/chimera/counts_matrix.tsv.gz",
        "results/chimera/candidates.tsv.gz",
    ]
    files += [
        f"results/chimera/read_evidence/per_sample/{s}_junction_qc.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/qc/junction_qc_mqc.json",
        "results/chimera/qc/te_gene_chimeras_mqc.json",
        "results/chimera/qc/canonical_rate_mqc.json",
        "results/chimera/qc/junction_highlights_mqc.json",
        "results/chimera/qc/chimera_candidates_table_mqc.json",
        "results/chimera/qc/chimera_evidence_guide_mqc.json",
        "results/chimera/qc/chimera_evidence_composition_mqc.json",
        "results/chimera/qc/chimera_evidence_correlation_mqc.json",
        "results/chimera/qc/chimera_evidence_candidates_mqc.json",
        "results/chimera/qc/sample_evidence_status_mqc.json",
    ]
    if WRITE_IGV_BED:
        files += [
            f"results/chimera/read_evidence/igv/{s}_junctions.bed"
            for s in SAMPLES
        ]
    return files


def all_sample_qc_outputs():
    """Sample-QC artifacts for the `all` target (Snakefile). Only produced
    when a counts matrix is written (it is the QC view's input)."""
    if not WRITE_COUNTS:
        return []
    transform = CHIMERA_QC["pca_transform"]
    return [
        f"results/chimera/qc/{transform}_counts.tsv.gz",
        f"results/chimera/qc/pca_{transform}_mqc.json",
        f"results/chimera/qc/heatmap_{transform}_mqc.json",
    ]


rule parse_chimeric_junctions:
    # Annotates one sample's STAR chimeric junctions against the gene/TE
    # tracks. See parse_chimeric_junctions.py for the full column spec.
    # Pure-python, so it runs in the base environment.
    input:
        junctions="results/star/{sample}_Chimeric.out.junction",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
        strandedness=strandedness_input,
    output:
        junctions="results/chimera/read_evidence/per_sample/{sample}_junctions.tsv.gz",
        te_gene_chimeras="results/chimera/read_evidence/per_sample/{sample}_junctions_te-gene-chimeras.tsv.gz",
    params:
        tolerance=config["chimera"]["junction"]["breakpoint_tolerance"],
        canonical_flag=lambda wc: (
            "--require-canonical"
            if config["chimera"]["junction"]["require_canonical_junction"]
            else ""
        ),
        library=get_strandedness_param,
    threads: get_resources("parse_chimeric_junctions")["threads"]
    resources:
        mem_mb=get_resources("parse_chimeric_junctions")["mem_mb"],
        runtime=get_resources("parse_chimeric_junctions")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/parse_chimeric_junctions/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/parse/{sample}.log",
    shell:
        "python3 {SCRIPTS_DIR}/parse_chimeric_junctions.py "
        "--junctions {input.junctions} "
        "--genes {input.genes} --exons {input.exons} --te {input.te} "
        "--sample {wildcards.sample} "
        "--breakpoint-tolerance {params.tolerance} "
        "{params.canonical_flag} "
        "--library-strandedness {params.library} "
        "--out {output.junctions} "
        "--te-out {output.te_gene_chimeras} > {log} 2>&1"


def _telocal_counts_for_chimera(wc=None):
    """Per-sample TElocal cntTable paths (only when TELOCAL_ENABLED)."""
    return [f"results/telocal/{s}.cntTable.gz" for s in SAMPLES]


rule chimera_telocal_index:
    # Builds the merged, all-samples TElocal interval index ONCE (sparse loci
    # won't overlap if a sample lacks them, which is fine) instead of every
    # per-sample chimera_telocal_annotate job re-parsing and re-summing all
    # samples' cntTables from scratch -- an Nx redundant rebuild of the
    # identical index for N samples.  Mirrors telocal_locind's
    # build-once-reuse-everywhere shape.  See chimera_telocal_index.py /
    # build_chimera_telocal_index.py for the compact columnar representation.
    input:
        telocal_tables=_telocal_counts_for_chimera,
        # The coordinate source. A TElocal cntTable key carries coordinates
        # only when the TE GTF's transcript_id happens to be a coordinate
        # string; otherwise this BED is the only way to place a locus.
        locations="results/telocal/telocal_locations.bed",
    output:
        "results/chimera/read_evidence/telocal_index.pkl.gz",
    threads: get_resources("chimera_telocal_index")["threads"]
    resources:
        mem_mb=get_resources("chimera_telocal_index")["mem_mb"],
        runtime=get_resources("chimera_telocal_index")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_telocal_index/chimera_telocal_index.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/telocal_index.log",
    shell:
        "python3 {SCRIPTS_DIR}/build_chimera_telocal_index.py "
        "--telocal-tables {input.telocal_tables} "
        "--locations {input.locations} "
        "--out {output} > {log} 2>&1"


rule chimera_telocal_annotate:
    # Enrich per-sample junction tables with TElocal locus-level expression.
    # Adds telocal_count, telocal_locus, telocal_active, te_refined_by_telocal
    # columns.  Only runs when both chimera and telocal are enabled.
    # Consumes the shared index built once by chimera_telocal_index -- this
    # job's own memory need no longer scales with sample count (it's the
    # same fixed index regardless of cohort size), unlike the old
    # per-sample-rebuild version.
    input:
        junctions="results/chimera/read_evidence/per_sample/{sample}_junctions.tsv.gz",
        telocal_index="results/chimera/read_evidence/telocal_index.pkl.gz",
    output:
        junctions="results/chimera/read_evidence/per_sample/{sample}_junctions_with-telocal.tsv.gz",
        te_gene_chimeras=(
            "results/chimera/read_evidence/per_sample/{sample}_junctions_te-gene-chimeras_with-telocal.tsv.gz"
        ),
    params:
        sample_name=lambda wc: wc.sample,
        tolerance=config["chimera"]["junction"]["breakpoint_tolerance"],
    threads: get_resources("chimera_telocal_annotate")["threads"]
    resources:
        mem_mb=get_resources("chimera_telocal_annotate")["mem_mb"],
        runtime=get_resources("chimera_telocal_annotate")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_telocal_annotate/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/telocal_annotate/{sample}.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_telocal_annotate.py "
        "--junctions {input.junctions} "
        "--telocal-index {input.telocal_index} "
        "--sample-name {params.sample_name} "
        "--breakpoint-tolerance {params.tolerance} "
        "--out {output.junctions} "
        "--te-out {output.te_gene_chimeras} > {log} 2>&1"


rule chimera_counts:
    # Merges every sample's junction table into the all-events catalog and
    # the event x sample counts matrix (chimera_counts.py). The counts matrix
    # feeds the sample-QC PCA/clustering stage.
    input:
        tables=chimera_counts_input(),
    output:
        events="results/chimera/read_evidence/all_events.tsv.gz",
        counts="results/chimera/counts_matrix.tsv.gz",
        te_events="results/chimera/read_evidence/te-gene-chimeras.tsv.gz",
    params:
        sample_names=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("chimera_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_counts"),
        runtime=get_resources("chimera_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_counts/chimera_counts.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/counts.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_counts.py "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--out-events {output.events} "
        "--out-counts {output.counts} "
        "--out-te-events {output.te_events} > {log} 2>&1"


rule junction_qc:
    # Per-sample junction QC summary (junction_qc.py) for the MultiQC custom
    # content table.
    input:
        "results/chimera/read_evidence/per_sample/{sample}_junctions.tsv.gz",
    output:
        "results/chimera/read_evidence/per_sample/{sample}_junction_qc.tsv.gz",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("junction_qc")["threads"]
    resources:
        mem_mb=get_resources("junction_qc")["mem_mb"],
        runtime=get_resources("junction_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/junction_qc/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/junction_qc/{sample}.log",
    shell:
        "python3 {SCRIPTS_DIR}/junction_qc.py "
        "--table {input} --sample {params.sample} --out {output} > {log} 2>&1"


rule chimera_evidence:
    # The unified gene-TE candidate catalogue: one row per (gene, TE
    # insertion) pair with every line of evidence either screen produced, and
    # no score (chimera_evidence.py -- the confidence tier it used to carry
    # was removed for lacking a validated weighting).  The two screens key
    # their output differently
    # -- junction rows are breakpoint-keyed, assembly rows transcript-keyed
    # -- so this is the only place they can be read side by side.
    #
    # Lives in the junction rule file (rather than the assembly one) because
    # the junction screen is the anchor: it always runs when chimera
    # detection is on, while the assembly screen is a separate switch.  Its
    # candidates.tsv.gz is therefore an OPTIONAL input, and the table simply
    # reports found_by: junction for everything when it is absent.
    input:
        junction="results/chimera/read_evidence/te-gene-chimeras.tsv.gz",
        **({"assembly": "results/chimera/transcript_evidence/transcripts.tsv.gz"}
           if CHIMERA_ASSEMBLY_ENABLED else {}),
    output:
        "results/chimera/candidates.tsv.gz",
    params:
        assembly=(
            "--assembly results/chimera/transcript_evidence/transcripts.tsv.gz"
            if CHIMERA_ASSEMBLY_ENABLED else ""
        ),
    threads: get_resources("chimera_evidence")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_evidence"),
        runtime=get_resources("chimera_evidence")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_evidence/chimera_evidence.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/chimera_evidence.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_evidence.py "
        "--junction {input.junction} {params.assembly} "
        "--out {output} > {log} 2>&1"


rule sample_evidence_status:
    # FastQC-style status grid: samples down the side, evidence layers
    # across the top (sample_evidence_status_mqc.py).
    #
    # The other two evidence heatmaps are cohort-level and cannot answer the
    # first question after a run: is ONE of my samples different? A replicate
    # contributing a tenth of the gene-TE junctions its peers do is absorbed
    # by any cohort summary. Scored against the cohort median rather than
    # absolute cut-offs, since chimeric yield varies by orders of magnitude
    # with depth and library prep.
    input:
        qc_tables=lambda wc: [
            f"results/chimera/read_evidence/per_sample/{s}_junction_qc.tsv.gz" for s in SAMPLES
        ],
        strandedness="results/rseqc/strandedness_check_mqc.json",
    output:
        "results/chimera/qc/sample_evidence_status_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("sample_evidence_status")["threads"]
    resources:
        mem_mb=get_resources("sample_evidence_status")["mem_mb"],
        runtime=get_resources("sample_evidence_status")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/sample_evidence_status/sample_evidence_status.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/sample_evidence_status.log",
    shell:
        "python3 {SCRIPTS_DIR}/sample_evidence_status_mqc.py "
        "--qc-tables {input.qc_tables} --samples {params.samples} "
        "--strandedness {input.strandedness} "
        "--out {output} > {log} 2>&1"


rule chimera_candidates_table:
    # The report's candidate list. A table, not a ranking: MultiQC's native
    # table sorts on any column, so the reader supplies the ordering the
    # pipeline deliberately does not
    # (chimera_candidates_table_mqc.py explains why).
    input:
        evidence="results/chimera/candidates.tsv.gz",
        gene_names="results/reference/gene_id_to_name.tsv.gz",
    output:
        "results/chimera/qc/chimera_candidates_table_mqc.json",
    params:
        top_n=CHIMERA_TABLE_TOP_N,
        source_path="results/chimera/candidates.tsv.gz",
    threads: get_resources("chimera_candidates_table")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_candidates_table"),
        runtime=get_resources("chimera_candidates_table")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_candidates_table/chimera_candidates_table.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/candidates_table.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_candidates_table_mqc.py "
        "--evidence {input.evidence} --gene-names {input.gene_names} "
        "--top-n {params.top_n} --source-path {params.source_path} "
        "--out {output} > {log} 2>&1"


rule chimera_evidence_guide:
    # The report's guide to reading the chimera evidence -- what each signal
    # is worth and what has been measured about it -- plus the cohort's
    # evidence composition (chimera_evidence_guide_mqc.py).
    #
    # Deliberately renders NO candidate table. This section replaced a
    # four-tier confidence ladder that was removed for lacking any validated
    # weighting; rendering even an unranked top-N would re-create it, because
    # whatever order the rows land in reads as importance. The full catalogue
    # is candidates.tsv.gz, which the reader sorts for their own question.
    #
    # No gene_id -> gene_name input: with no table there is nothing to label.
    input:
        evidence="results/chimera/candidates.tsv.gz",
    output:
        guide="results/chimera/qc/chimera_evidence_guide_mqc.json",
        composition="results/chimera/qc/chimera_evidence_composition_mqc.json",
    threads: get_resources("chimera_evidence_guide")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_evidence_guide"),
        runtime=get_resources("chimera_evidence_guide")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_evidence_guide/chimera_evidence_guide.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/evidence_guide.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_evidence_guide_mqc.py "
        "--evidence {input.evidence} "
        "--out-guide {output.guide} "
        "--out-composition {output.composition} > {log} 2>&1"


rule chimera_evidence_heatmap:
    # Two heatmaps over chimera_evidence.tsv.gz and NO score
    # (chimera_evidence_heatmap.py): dimension x dimension correlation, and
    # the union of the per-dimension leaders as a candidate view.
    #
    # Exists because collapsing these measurements into one ordinal tier hid
    # the opposite of what was assumed, twice -- TElocal expression is
    # anti-correlated with the splice motif, and screen agreement sits near
    # its chance rate. Both were plain the moment the dimensions were shown
    # against each other. Judge the evidence, then decide on a ranking.
    input:
        evidence="results/chimera/candidates.tsv.gz",
        gene_names="results/reference/gene_id_to_name.tsv.gz",
    output:
        correlation="results/chimera/qc/chimera_evidence_correlation_mqc.json",
        candidates="results/chimera/qc/chimera_evidence_candidates_mqc.json",
    threads: get_resources("chimera_evidence_heatmap")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_evidence_heatmap"),
        runtime=get_resources("chimera_evidence_heatmap")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_evidence_heatmap/chimera_evidence_heatmap.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/evidence_heatmap.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_evidence_heatmap.py "
        "--evidence {input.evidence} --gene-names {input.gene_names} "
        "--out-correlation {output.correlation} "
        "--out-candidates {output.candidates} > {log} 2>&1"


rule junction_highlights:
    # The junction screen's "what to look at first" guide + ranked top-N
    # table (junction_highlights_mqc.py) -- the mirror of the assembly
    # screen's highlights section, so a default run's report has a guided
    # entry point for BOTH screens rather than only the assembly one.
    # Reads the merged gene<->TE table (not the per-sample QC metrics),
    # because the ranking needs per-event annotation columns.
    input:
        te_events="results/chimera/read_evidence/te-gene-chimeras.tsv.gz",
    output:
        "results/chimera/qc/junction_highlights_mqc.json",
    threads: get_resources("junction_highlights")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("junction_highlights"),
        runtime=get_resources("junction_highlights")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/junction_highlights/junction_highlights.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/junction_highlights.log",
    shell:
        "python3 {SCRIPTS_DIR}/junction_highlights_mqc.py "
        "--te-events {input.te_events} "
        "--out {output} > {log} 2>&1"


rule junction_qc_barplot:
    # Merges the per-sample junction QC tables into two MultiQC bar-plot
    # custom-content documents (junction_qc_mqc.py): per-sample direction
    # composition as counts and % of total junctions, plus the gene<->TE
    # subset (the gene-TE chimeras view), rendered inside multiqc_report.html in
    # the custom_content module.
    input:
        tables=lambda wc: [
            f"results/chimera/read_evidence/per_sample/{s}_junction_qc.tsv.gz" for s in SAMPLES
        ],
    output:
        junction="results/chimera/qc/junction_qc_mqc.json",
        te_gene_chimeras="results/chimera/qc/te_gene_chimeras_mqc.json",
        canonical="results/chimera/qc/canonical_rate_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("junction_qc_barplot")["threads"]
    resources:
        mem_mb=get_resources("junction_qc_barplot")["mem_mb"],
        runtime=get_resources("junction_qc_barplot")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/junction_qc_barplot/junction_qc_barplot.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/junction_qc_barplot.log",
    shell:
        "python3 {SCRIPTS_DIR}/junction_qc_mqc.py "
        "--tables {input.tables} "
        "--samples {params.samples} "
        "--out {output.junction} "
        "--out-te-gene-chimeras {output.te_gene_chimeras} "
        "--out-canonical {output.canonical} > {log} 2>&1"


rule chimera_igv_bed:
    # Optional IGV track per sample (one BED12-style row per junction), so
    # candidates can be inspected visually. Gated by
    # config["chimera"]["junction"]["outputs"]["write_igv_bed"].
    input:
        "results/chimera/read_evidence/per_sample/{sample}_junctions.tsv.gz",
    output:
        "results/chimera/read_evidence/igv/{sample}_junctions.bed",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("chimera_igv_bed")["threads"]
    resources:
        mem_mb=get_resources("chimera_igv_bed")["mem_mb"],
        runtime=get_resources("chimera_igv_bed")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_igv_bed/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_junction/igv/{sample}.log",
    script:
        "../scripts/junctions_to_igv_bed.py"
