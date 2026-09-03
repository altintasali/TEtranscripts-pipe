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
# counts matrix the sample-QC stage (chimera_reads_qc.smk) consumes.
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
#   chimera_reads_classify per-sample junction annotation (+ the gene-TE subset)
#   chimera_reads_counts           merge per-sample tables -> all_events + counts
#                            + te-gene-chimeras
#   chimera_reads_qc              per-sample QC summary (MultiQC custom content)
#   chimera_reads_igv_bed          per-sample IGV track (config-gated)
# -----------------------------------------------------------------------------
import os

WRITE_COUNTS = bool(config["chimera"]["reads"]["outputs"]["write_counts_matrix"])
WRITE_IGV_BED = bool(config["chimera"]["reads"]["outputs"]["write_igv_bed"])
# CHIMERA_QC now lives in common/runtime.smk -- the assembly screen needs it
# too, and this file is not included when the junction screen is off.


def chimera_reads_counts_input():
    if TELOCAL_ENABLED:
        return [
            f"results/chimera/reads/per_sample/{s}_junctions_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
    return [
        f"results/chimera/reads/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]


def all_chimera_outputs():
    """Chimera artifacts for the `all` target (Snakefile)."""
    files = [
        f"results/chimera/reads/per_sample/{s}_junctions.tsv.gz"
        for s in SAMPLES
    ]
    if TELOCAL_ENABLED:
        files += [
            f"results/chimera/reads/per_sample/{s}_junctions_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
        files += [
            f"results/chimera/reads/per_sample/{s}_junctions_te-gene-chimeras_with-telocal.tsv.gz"
            for s in SAMPLES
        ]
    files += [
        f"results/chimera/reads/per_sample/{s}_junctions_te-gene-chimeras.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/reads/all_events.tsv.gz",
        "results/chimera/reads/te-gene-chimeras.tsv.gz",
        "results/chimera/counts_matrix.tsv.gz",
        "results/chimera/candidates.tsv.gz",
        "results/chimera/candidates_explorer.html",
    ]
    files += [
        f"results/chimera/reads/per_sample/{s}_chimera_reads_qc.tsv.gz"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/qc/chimera_reads_qc_mqc.json",
        "results/chimera/qc/te_gene_chimeras_mqc.json",
        "results/chimera/qc/canonical_rate_mqc.json",
        "results/chimera/qc/chimera_reads_highlights_mqc.json",
        "results/chimera/qc/chimera_candidates_table_mqc.json",
        "results/chimera/qc/chimera_reads_te_type_mqc.json",
        "results/chimera/qc/chimera_evidence_guide_mqc.json",
        "results/chimera/qc/chimera_evidence_composition_mqc.json",
        "results/chimera/qc/chimera_evidence_correlation_mqc.json",
        "results/chimera/qc/chimera_evidence_candidates_mqc.json",
    ]
    if WRITE_IGV_BED:
        files += [
            f"results/chimera/reads/igv/{s}_junctions.bed"
            for s in SAMPLES
        ]
    return files


def all_chimera_reads_sample_qc_outputs():
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


rule chimera_reads_classify:
    # Annotates one sample's STAR chimeric junctions against the gene/TE
    # tracks. See classify_chimera_reads.py for the full column spec.
    # Pure-python, so it runs in the base environment.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/classify_chimera_reads.py",
        junctions="results/star/{sample}_Chimeric.out.junction",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
        strandedness=strandedness_input,
    output:
        junctions="results/chimera/reads/per_sample/{sample}_junctions.tsv.gz",
        te_gene_chimeras="results/chimera/reads/per_sample/{sample}_junctions_te-gene-chimeras.tsv.gz",
    params:
        tolerance=config["chimera"]["reads"]["breakpoint_tolerance"],
        canonical_flag=lambda wc: (
            "--require-canonical"
            if config["chimera"]["reads"]["require_canonical_junction"]
            else ""
        ),
        library=get_strandedness_param,
    threads: get_resources("chimera_reads_classify")["threads"]
    resources:
        mem_mb=get_resources("chimera_reads_classify")["mem_mb"],
        runtime=get_resources("chimera_reads_classify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_classify/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/classify/{sample}.log",
    shell:
        "python3 {input.script} "
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
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/build_chimera_telocal_index.py",
        telocal_tables=_telocal_counts_for_chimera,
        # The coordinate source. A TElocal cntTable key carries coordinates
        # only when the TE GTF's transcript_id happens to be a coordinate
        # string; otherwise this BED is the only way to place a locus.
        locations="results/telocal/telocal_locations.bed",
    output:
        "results/chimera/reads/telocal_index.pkl.gz",
    threads: get_resources("chimera_telocal_index")["threads"]
    resources:
        mem_mb=get_resources("chimera_telocal_index")["mem_mb"],
        runtime=get_resources("chimera_telocal_index")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_telocal_index/chimera_telocal_index.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/telocal_index.log",
    shell:
        "python3 {input.script} "
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
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_telocal_annotate.py",
        junctions="results/chimera/reads/per_sample/{sample}_junctions.tsv.gz",
        telocal_index="results/chimera/reads/telocal_index.pkl.gz",
    output:
        junctions="results/chimera/reads/per_sample/{sample}_junctions_with-telocal.tsv.gz",
        te_gene_chimeras=(
            "results/chimera/reads/per_sample/{sample}_junctions_te-gene-chimeras_with-telocal.tsv.gz"
        ),
    params:
        sample_name=lambda wc: wc.sample,
        tolerance=config["chimera"]["reads"]["breakpoint_tolerance"],
    threads: get_resources("chimera_telocal_annotate")["threads"]
    resources:
        mem_mb=get_resources("chimera_telocal_annotate")["mem_mb"],
        runtime=get_resources("chimera_telocal_annotate")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_telocal_annotate/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/telocal_annotate/{sample}.log",
    shell:
        "python3 {input.script} "
        "--junctions {input.junctions} "
        "--telocal-index {input.telocal_index} "
        "--sample-name {params.sample_name} "
        "--breakpoint-tolerance {params.tolerance} "
        "--out {output.junctions} "
        "--te-out {output.te_gene_chimeras} > {log} 2>&1"


rule chimera_reads_counts:
    # Merges every sample's junction table into the all-events catalog and
    # the event x sample counts matrix (chimera_reads_counts.py). The counts matrix
    # feeds the sample-QC PCA/clustering stage.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_reads_counts.py",
        tables=chimera_reads_counts_input(),
    output:
        events="results/chimera/reads/all_events.tsv.gz",
        counts="results/chimera/counts_matrix.tsv.gz",
        te_events="results/chimera/reads/te-gene-chimeras.tsv.gz",
    params:
        sample_names=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("chimera_reads_counts")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_reads_counts"),
        runtime=get_resources("chimera_reads_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_counts/chimera_reads_counts.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/counts.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--out-events {output.events} "
        "--out-counts {output.counts} "
        "--out-te-events {output.te_events} > {log} 2>&1"


rule chimera_reads_qc:
    # Per-sample junction QC summary (chimera_reads_qc.py) for the MultiQC custom
    # content table.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_reads_qc.py",
        table="results/chimera/reads/per_sample/{sample}_junctions.tsv.gz",
    output:
        "results/chimera/reads/per_sample/{sample}_chimera_reads_qc.tsv.gz",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("chimera_reads_qc")["threads"]
    resources:
        mem_mb=get_resources("chimera_reads_qc")["mem_mb"],
        runtime=get_resources("chimera_reads_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_qc/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_reads_qc/{sample}.log",
    shell:
        "python3 {input.script} "
        "--table {input.table} --sample {params.sample} --out {output} > {log} 2>&1"


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
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_evidence.py",
        junction="results/chimera/reads/te-gene-chimeras.tsv.gz",
        **({"assembly": "results/chimera/assembly/transcripts.tsv.gz"}
           if CHIMERA_ASSEMBLY_ENABLED else {}),
    output:
        "results/chimera/candidates.tsv.gz",
    params:
        assembly=(
            "--assembly results/chimera/assembly/transcripts.tsv.gz"
            if CHIMERA_ASSEMBLY_ENABLED else ""
        ),
    threads: get_resources("chimera_evidence")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_evidence"),
        runtime=get_resources("chimera_evidence")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_evidence/chimera_evidence.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_evidence.log",
    shell:
        "python3 {input.script} "
        "--junction {input.junction} {params.assembly} "
        "--out {output} > {log} 2>&1"


rule chimera_reads_te_type:
    # The read screen's TE-type view, per sample. Its counterpart is the
    # assembly screen's own class chart (Assembly - TE type); the two are kept
    # SEPARATE because they use the same three words for different
    # measurements -- position relative to the gene here, transcript exon
    # structure there -- and one shared plot invites reading agreement as
    # corroboration. Each section says so.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_reads_te_type_mqc.py",
        qc_tables=expand("results/chimera/reads/per_sample/"
                         "{sample}_chimera_reads_qc.tsv.gz", sample=SAMPLES),
    output:
        "results/chimera/qc/chimera_reads_te_type_mqc.json",
    params:
        samples=SAMPLES,
    threads: get_resources("chimera_reads_te_type")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_reads_te_type"),
        runtime=get_resources("chimera_reads_te_type")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_te_type/chimera_reads_te_type.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/te_type.log",
    shell:
        "python3 {input.script} "
        "--qc-tables {input.qc_tables} --samples {params.samples} "
        "--out {output} > {log} 2>&1"


rule chimera_candidates_table:
    # The report's candidate list. A table, not a ranking: MultiQC's native
    # table sorts on any column, so the reader supplies the ordering the
    # pipeline deliberately does not
    # (chimera_candidates_table_mqc.py explains why).
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_candidates_table_mqc.py",
        evidence="results/chimera/candidates.tsv.gz",
        gene_names="results/reference/gene_id_to_name.tsv.gz",
    output:
        "results/chimera/qc/chimera_candidates_table_mqc.json",
    params:
        top_n=CHIMERA_TABLE_TOP_N,
        source_path="results/chimera/candidates.tsv.gz",
        explorer_path="results/chimera/candidates_explorer.html",
    threads: get_resources("chimera_candidates_table")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_candidates_table"),
        runtime=get_resources("chimera_candidates_table")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_candidates_table/chimera_candidates_table.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/candidates_table.log",
    shell:
        "python3 {input.script} "
        "--evidence {input.evidence} --gene-names {input.gene_names} "
        "--top-n {params.top_n} --source-path {params.source_path} "
        "--explorer-path {params.explorer_path} "
        "--out {output} > {log} 2>&1"


rule chimera_candidates_explorer:
    # Standalone, self-contained interactive HTML over the FULL candidates
    # catalogue -- the top_n cap on chimera_candidates_table above exists
    # because MultiQC embeds table data into one already-large report; this
    # is the "all rows, sortable/searchable/filterable, no MultiQC, no
    # server" companion (chimera_candidates_explorer.R, DT + htmlwidgets).
    # Also adds two columns candidates.tsv.gz does not carry: a ready-to-
    # paste IGV locus per gene and per TE, joined from genes.bed/te.bed
    # (both keyed on the same gene_id/te_id candidates.tsv.gz uses, and
    # always written whenever a chimera screen is enabled -- unlike the
    # optional, differently-keyed chimera_reads_igv_bed/
    # chimera_assembly_igv_bed tracks, which cannot be searched by
    # candidate name).
    #
    # genes.bed/te.bed are GENOME-WIDE -- every gene and every TE insertion
    # in the annotation, not just the ones with a candidate. te.bed
    # especially: a real mouse/human TE annotation is millions of rows.
    # Loading that whole file into R OOM-killed this rule on real data (2.3
    # GB against a 2.4 GB request, measured on a real run -- a synthetic
    # dev fixture sized to just the candidate count never caught this).
    # awk-filtering both BEDs down to only the gene_id/te_id values
    # candidates.tsv.gz actually references, BEFORE R ever sees them, is
    # the same fix as rseqc_gene_body_coverage's BED thinning in
    # bam_qc.smk: cheap, streaming, and keeps R's peak memory proportional
    # to candidate count instead of genome size.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_candidates_explorer.R",
        evidence="results/chimera/candidates.tsv.gz",
        gene_names="results/reference/gene_id_to_name.tsv.gz",
        genes_bed="results/reference/genes.bed",
        te_bed="results/reference/te.bed",
    output:
        "results/chimera/candidates_explorer.html",
    threads: get_resources("chimera_candidates_explorer")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_candidates_explorer"),
        runtime=get_resources("chimera_candidates_explorer")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_candidates_explorer/chimera_candidates_explorer.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/candidates_explorer.log",
    conda:
        CANDIDATES_EXPLORER_ENV
    shell:
        "genes_ids={resources.tmpdir}/candidates_gene_ids.txt; "
        "te_ids={resources.tmpdir}/candidates_te_ids.txt; "
        "genes_filtered={resources.tmpdir}/candidates_genes.bed; "
        "te_filtered={resources.tmpdir}/candidates_te.bed; "
        "gzip -dc {input.evidence} | tail -n +2 | cut -f1 | sort -u > \"$genes_ids\" && "
        "gzip -dc {input.evidence} | tail -n +2 | cut -f2 | sort -u > \"$te_ids\" && "
        "awk -F'\\t' 'NR==FNR{{ids[$1]=1; next}} ($4 in ids)' "
        "\"$genes_ids\" {input.genes_bed} > \"$genes_filtered\" && "
        "awk -F'\\t' 'NR==FNR{{ids[$1]=1; next}} ($4 in ids)' "
        "\"$te_ids\" {input.te_bed} > \"$te_filtered\" && "
        "Rscript {input.script} "
        "{input.evidence} {input.gene_names} \"$genes_filtered\" \"$te_filtered\" "
        "{output} > {log} 2>&1"


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
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_evidence_guide_mqc.py",
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
        "results/pipeline_info/logs/chimera_reads/evidence_guide.log",
    shell:
        "python3 {input.script} "
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
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_evidence_heatmap.py",
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
        "results/pipeline_info/logs/chimera_reads/evidence_heatmap.log",
    shell:
        "python3 {input.script} "
        "--evidence {input.evidence} --gene-names {input.gene_names} "
        "--out-correlation {output.correlation} "
        "--out-candidates {output.candidates} > {log} 2>&1"


rule chimera_reads_highlights:
    # The junction screen's "what to look at first" guide + ranked top-N
    # table (chimera_reads_highlights_mqc.py) -- the mirror of the assembly
    # screen's highlights section, so a default run's report has a guided
    # entry point for BOTH screens rather than only the assembly one.
    # Reads the merged gene<->TE table (not the per-sample QC metrics),
    # because the ranking needs per-event annotation columns.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_reads_highlights_mqc.py",
        te_events="results/chimera/reads/te-gene-chimeras.tsv.gz",
    output:
        "results/chimera/qc/chimera_reads_highlights_mqc.json",
    threads: get_resources("chimera_reads_highlights")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_reads_highlights"),
        runtime=get_resources("chimera_reads_highlights")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_highlights/chimera_reads_highlights.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_reads_highlights.log",
    shell:
        "python3 {input.script} "
        "--te-events {input.te_events} "
        "--out {output} > {log} 2>&1"


rule chimera_reads_qc_barplot:
    # Merges the per-sample junction QC tables into two MultiQC bar-plot
    # custom-content documents (chimera_reads_qc_mqc.py): per-sample direction
    # composition as counts and % of total junctions, plus the gene<->TE
    # subset (the gene-TE chimeras view), rendered inside multiqc_report.html in
    # the custom_content module.
    input:
        # Declared so that EDITING the script re-runs the rule.
        # Snakemake's code trigger hashes the shell command STRING,
        # not the file it names, so without this an edit to the
        # script leaves stale outputs in place silently.
        script=f"{SCRIPTS_DIR}/chimera_reads_qc_mqc.py",
        tables=lambda wc: [
            f"results/chimera/reads/per_sample/{s}_chimera_reads_qc.tsv.gz" for s in SAMPLES
        ],
    output:
        junction="results/chimera/qc/chimera_reads_qc_mqc.json",
        te_gene_chimeras="results/chimera/qc/te_gene_chimeras_mqc.json",
        canonical="results/chimera/qc/canonical_rate_mqc.json",
        enrichment="results/chimera/qc/canonical_enrichment_mqc.json",
    params:
        samples=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("chimera_reads_qc_barplot")["threads"]
    resources:
        mem_mb=get_resources("chimera_reads_qc_barplot")["mem_mb"],
        runtime=get_resources("chimera_reads_qc_barplot")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_qc_barplot/chimera_reads_qc_barplot.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/chimera_reads_qc_barplot.log",
    shell:
        "python3 {input.script} "
        "--tables {input.tables} "
        "--samples {params.samples} "
        "--out {output.junction} "
        "--out-te-gene-chimeras {output.te_gene_chimeras} "
        "--out-canonical {output.canonical} "
        "--out-enrichment {output.enrichment} > {log} 2>&1"


rule chimera_reads_igv_bed:
    # Optional IGV track per sample (one BED12-style row per junction), so
    # candidates can be inspected visually. Gated by
    # config["chimera"]["reads"]["outputs"]["write_igv_bed"].
    input:
        "results/chimera/reads/per_sample/{sample}_junctions.tsv.gz",
    output:
        "results/chimera/reads/igv/{sample}_junctions.bed",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("chimera_reads_igv_bed")["threads"]
    resources:
        mem_mb=get_resources("chimera_reads_igv_bed")["mem_mb"],
        runtime=get_resources("chimera_reads_igv_bed")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_reads_igv_bed/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_reads/igv/{sample}.log",
    script:
        "../scripts/chimera_reads_to_igv_bed.py"
