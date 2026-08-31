# -----------------------------------------------------------------------------
# Chimera-assembly screen: gene-TE chimera detection from StringTie assembly
# structure, complementing chimera_junction.smk's read-level screen.
#
# EXPERIMENTAL and OFF BY DEFAULT: newer and less validated than
# chimera_junction. STAR only flags a junction as "chimeric" when a read
# can't be explained by one linear (possibly spliced) alignment -- a TE that
# splices into a gene via an ordinary, canonical, nearby intron aligns as a
# completely normal spliced read and never reaches chimera_junction at all.
# This screen catches that case instead, from StringTie's assembled
# transcript structure.
#
# Rules:
#   star_align_for_assembly  dedicated per-sample STAR pass (--outSAMstrandField
#                             intronMotif, for StringTie's unstranded-data
#                             strand inference) -- feeds ONLY this screen; the
#                             main alignment everything else uses is untouched.
#   stringtie_assemble        per-sample de novo assembly
#   stringtie_merge           cross-sample structural union
#   stringtie_requantify       per-sample re-quantification (-e -B) for TPM
#   chimera_assembly_classify  structural classification (te_initiated/
#                               te_exonized/te_terminated/unspliced_te_only)
#   chimera_assembly_quantify  candidate x sample TPM matrix
#   chimera_assembly_cross_evidence  cross-check against chimera_junction's calls
#
# genes.bed/exons.bed/te.bed are built by ref.smk's annotation_to_bed rule
# (shared with chimera_junction.smk).
# -----------------------------------------------------------------------------
import os

WRITE_IGV_BED_ASSEMBLY = bool(config["chimera"]["assembly"]["outputs"]["write_igv_bed"])


def all_chimera_assembly_outputs():
    """Chimera-assembly artifacts for the `all` target (Snakefile)."""
    files = [
        "results/chimera/transcript_evidence/transcripts.tsv.gz",
        "results/chimera/transcript_evidence/tpm_matrix.tsv.gz",
        "results/chimera/qc/chimera_assembly_classes_mqc.json",
        "results/chimera/qc/chimera_assembly_highlights_mqc.json",
        "results/chimera/qc/chimera_assembly_strand_rate_mqc.json",
        # the PCA/Clusters view, matching the read screen's
        "results/chimera/qc/assembly_pca_log2_mqc.json",
        "results/chimera/qc/assembly_heatmap_log2_mqc.json",
    ]
    if CHIMERA_JUNCTION_ENABLED:
        files.append("results/chimera/transcript_evidence/transcripts_with_read_support.tsv.gz")
    if WRITE_IGV_BED_ASSEMBLY:
        files.append("results/chimera/transcript_evidence/igv/transcripts.bed")
    return files


def _chimera_assembly_summary_input():
    """Prefer the cross-referenced candidates table (adds
    confirmed_by_junction_screen) when the junction screen also ran; fall
    back to the plain candidates table otherwise."""
    if CHIMERA_JUNCTION_ENABLED:
        return "results/chimera/transcript_evidence/transcripts_with_read_support.tsv.gz"
    return "results/chimera/transcript_evidence/transcripts.tsv.gz"


def _assembly_two_pass_args(wildcards, input):
    """Same STAR 2-pass handoff as align.smk's star_align (config
    star.two_pass), so this dedicated pass isn't working from worse splice
    detection than the main alignment."""
    if STAR_TWO_PASS == "per_sample":
        return "--twopassMode Basic"
    if STAR_TWO_PASS == "cohort":
        return f"--sjdbFileChrStartEnd {input.merged_sj}"
    return ""


def stringtie_strand_flag(wildcards, input):
    """--fr / --rf / nothing, from the same per-sample strandedness
    resolution chimera_junction.smk/tetranscripts.smk already use.
    StringTie infers per-transcript strand for spliced reads from the XS
    tag (see star_align_for_assembly's --outSAMstrandField intronMotif)
    even on unstranded libraries -- so "no" still yields usable multi-exon
    transcript strand calls, just not for single-exon ones (those have no
    splice to carry an XS tag at all)."""
    lib = get_strandedness_param(wildcards, input)
    return {"forward": "--fr", "reverse": "--rf"}.get(lib, "")


rule star_align_for_assembly:
    # Dedicated, chimera_assembly-only STAR pass. --outSAMstrandField
    # intronMotif gives StringTie the XS strand tag it needs for
    # transcript-strand inference on unstranded data, which the main
    # alignment (used by the read-evidence screen/TEcount/TElocal) does not set.
    #
    # This is a real, bounded tradeoff, not a free win: per STAR's own
    # manual, that flag suppresses reads with non-canonical/undefined-
    # strand introns from its output -- so THIS PRIVATE BAM has reduced
    # sensitivity to non-canonical TE splices. The main BAM everything
    # else uses is completely untouched.
    input:
        unpack(star_input),
        idx=STAR_INDEX_DIR,
        **({"merged_sj": "results/star_pass1/merged_SJ.out.tab"}
           if STAR_TWO_PASS == "cohort" else {}),
    output:
        aln="results/chimera/transcript_evidence/per_sample/star/{sample}_Aligned.sortedByCoord.out.bam",
        log_final="results/chimera/transcript_evidence/per_sample/star/{sample}_Log.final.out",
    params:
        reads=star_reads_param,
        read_command=star_read_command_param,
        prefix=lambda wc: f"results/chimera/transcript_evidence/per_sample/star/{wc.sample}_",
        tmpdir=lambda wc: os.path.join(STAR_TMPDIR, f"star_assembly_{wc.sample}"),
        two_pass=_assembly_two_pass_args,
        extra=config["star"]["extra"],
    threads: get_resources("star_align_for_assembly")["threads"]
    resources:
        mem_mb=get_resources("star_align_for_assembly")["mem_mb"],
        runtime=get_resources("star_align_for_assembly")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/star_align_for_assembly/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/star/{sample}.log",
    conda:
        STAR_ENV
    shell:
        # Same benign-exit-crash tolerance as star_align (rules/align.smk).
        "mkdir -p results/chimera/transcript_evidence/per_sample/star && "
        "rm -rf {params.tmpdir} && "
        "trap 'rm -rf {params.tmpdir}' EXIT; "
        "(STAR --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {params.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outTmpDir {params.tmpdir} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMstrandField intronMotif "
        "{params.two_pass} "
        "{params.extra} "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether alignment output "
        "is actually complete anyway (see star_align rule comment re: known "
        "benign STAR exit-time crash)' >> {log}; "
        "if test -s {output.aln} && test -s {output.log_final} && "
        "grep -q 'ALL DONE!' {output.log_final}; then :; else "
        "echo '-- STAR failed for real; log tail --' >&2; "
        "tail -n 60 {log} >&2; exit 1; fi))"


rule stringtie_assemble:
    # Per-sample de novo assembly. -m 100 -c 1: loose thresholds (TEProf2's
    # recommendation) so low-coverage TE-driven transcripts aren't dropped
    # before classification -- filter on expression downstream instead.
    input:
        bam="results/chimera/transcript_evidence/per_sample/star/{sample}_Aligned.sortedByCoord.out.bam",
        strandedness=strandedness_input,
    output:
        gtf="results/chimera/transcript_evidence/per_sample/assembly/{sample}.gtf",
    params:
        strand=stringtie_strand_flag,
    threads: get_resources("stringtie_assemble")["threads"]
    resources:
        mem_mb=get_resources("stringtie_assemble")["mem_mb"],
        runtime=get_resources("stringtie_assemble")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/stringtie_assemble/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/assemble/{sample}.log",
    conda:
        STRINGTIE_ENV
    shell:
        "stringtie {input.bam} -o {output.gtf} "
        "-m 100 -c 1 {params.strand} -p {threads} > {log} 2>&1"


rule stringtie_merge:
    # Cross-sample structural union -- no -G, so genuinely novel TE-initiated
    # structures aren't discarded for lacking a matching annotated
    # transcript. Aggregates over every sample: the DAG sync point for this
    # screen (no sample's requantify starts until every sample's assembly
    # has finished), same shape as star_merge_junctions -- sized with
    # mem_per_sample from the start, not discovered the hard way.
    input:
        gtfs=expand("results/chimera/transcript_evidence/per_sample/assembly/{sample}.gtf", sample=SAMPLES),
    output:
        merged="results/chimera/transcript_evidence/stringtie_merge.gtf",
    threads: get_resources("stringtie_merge")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("stringtie_merge"),
        runtime=get_resources("stringtie_merge")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/stringtie_merge/stringtie_merge.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/merge.log",
    conda:
        STRINGTIE_ENV
    shell:
        "stringtie --merge -o {output.merged} {input.gtfs} > {log} 2>&1"


rule stringtie_requantify:
    # Per-sample re-quantification (-e -B) against the merged structure:
    # comparable per-sample TPM for every candidate transcript_id, the same
    # "assemble -> merge -> quantify" pattern nf-core/rnaseq's
    # --stringtie_ignore_gtf uses.
    input:
        bam="results/chimera/transcript_evidence/per_sample/star/{sample}_Aligned.sortedByCoord.out.bam",
        merged="results/chimera/transcript_evidence/stringtie_merge.gtf",
        strandedness=strandedness_input,
    output:
        gtf="results/chimera/transcript_evidence/per_sample/quant/{sample}.transcripts.gtf",
    params:
        strand=stringtie_strand_flag,
    threads: get_resources("stringtie_requantify")["threads"]
    resources:
        mem_mb=get_resources("stringtie_requantify")["mem_mb"],
        runtime=get_resources("stringtie_requantify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/stringtie_requantify/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/requantify/{sample}.log",
    conda:
        STRINGTIE_ENV
    shell:
        "stringtie -e -B -G {input.merged} -o {output.gtf} "
        "{params.strand} {input.bam} > {log} 2>&1"


rule chimera_assembly_classify:
    # Structural classification against the SAME genes.bed/exons.bed/te.bed
    # the junction screen uses (ref.smk's annotation_to_bed) -- no separate
    # reference-track build needed.
    input:
        gtf="results/chimera/transcript_evidence/stringtie_merge.gtf",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
    output:
        candidates="results/chimera/transcript_evidence/transcripts.tsv.gz",
    params:
        tolerance=config["chimera"]["assembly"]["breakpoint_tolerance"],
    threads: get_resources("chimera_assembly_classify")["threads"]
    resources:
        mem_mb=get_resources("chimera_assembly_classify")["mem_mb"],
        runtime=get_resources("chimera_assembly_classify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_assembly_classify/chimera_assembly_classify.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/classify.log",
    shell:
        "python3 {SCRIPTS_DIR}/classify_chimera_assembly.py "
        "--gtf {input.gtf} --genes {input.genes} --exons {input.exons} --te {input.te} "
        "--breakpoint-tolerance {params.tolerance} "
        "--out {output.candidates} > {log} 2>&1"


rule chimera_assembly_quantify:
    # Pulls per-sample TPM for every candidate transcript_id out of the
    # re-quantified GTFs -- an expression matrix to threshold/filter on
    # downstream (min TPM, min replicates), same "annotate-only, filter
    # later" philosophy as the junction screen.
    input:
        candidates="results/chimera/transcript_evidence/transcripts.tsv.gz",
        quant=expand("results/chimera/transcript_evidence/per_sample/quant/{sample}.transcripts.gtf", sample=SAMPLES),
    output:
        matrix="results/chimera/transcript_evidence/tpm_matrix.tsv.gz",
    params:
        sample_names=" ".join(SAMPLES),
    threads: get_resources("chimera_assembly_quantify")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_assembly_quantify"),
        runtime=get_resources("chimera_assembly_quantify")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_assembly_quantify/chimera_assembly_quantify.txt",
    log:
        "results/pipeline_info/logs/chimera_assembly/quantify.log",
    shell:
        "python3 {SCRIPTS_DIR}/quantify_chimera_assembly.py "
        "--candidates {input.candidates} --quant {input.quant} "
        "--sample-names {params.sample_names} --out {output.matrix} > {log} 2>&1"


if CHIMERA_JUNCTION_ENABLED:

    rule chimera_assembly_cross_evidence:
        # Cross-checks this screen's calls against the junction screen's
        # te-gene-chimeras.tsv.gz (STAR chimeric-junction-read-based): a
        # candidate backed by BOTH the read-level evidence AND the assembly
        # structure is much higher confidence than either alone. Only
        # meaningful (and only included) when the junction screen also runs.
        input:
            candidates="results/chimera/transcript_evidence/transcripts.tsv.gz",
            te_gene_chimeras="results/chimera/read_evidence/te-gene-chimeras.tsv.gz",
        output:
            "results/chimera/transcript_evidence/transcripts_with_read_support.tsv.gz",
        threads: get_resources("chimera_assembly_cross_evidence")["threads"]
        resources:
            mem_mb=get_resources("chimera_assembly_cross_evidence")["mem_mb"],
            runtime=get_resources("chimera_assembly_cross_evidence")["runtime"],
        log:
            "results/pipeline_info/logs/chimera_assembly/cross_evidence.log",
        shell:
            "python3 {SCRIPTS_DIR}/cross_evidence_chimera_assembly.py "
            "--candidates {input.candidates} --junction {input.te_gene_chimeras} "
            "--out {output} > {log} 2>&1"


rule chimera_assembly_summary_mqc:
    # MultiQC custom content: candidates-by-class barplot (split by
    # read-evidence confirmation when that screen ran too) plus this screen's
    # blind-spot note and the counts qualifying its own output. Uses the
    # cross-referenced table when available (see
    # _chimera_assembly_summary_input) so the confirmation split is possible.
    #
    # tpm_matrix is deliberately NOT an input: it was only ever read to sort
    # a ranked top-N table here, and the pipeline no longer ranks candidates
    # anywhere (see chimera_evidence_guide_mqc.py).
    input:
        candidates=_chimera_assembly_summary_input(),
    output:
        classes="results/chimera/qc/chimera_assembly_classes_mqc.json",
        highlights="results/chimera/qc/chimera_assembly_highlights_mqc.json",
        strand_rate="results/chimera/qc/chimera_assembly_strand_rate_mqc.json",
    threads: get_resources("chimera_assembly_summary_mqc")["threads"]
    resources:
        mem_mb=get_resources("chimera_assembly_summary_mqc")["mem_mb"],
        runtime=get_resources("chimera_assembly_summary_mqc")["runtime"],
    log:
        "results/pipeline_info/logs/chimera_assembly/summary_mqc.log",
    shell:
        "python3 {SCRIPTS_DIR}/chimera_assembly_summary_mqc.py "
        "--candidates {input.candidates} "
        "--out-classes {output.classes} --out-highlights {output.highlights} "
        "--out-strand-rate {output.strand_rate} "
        "> {log} 2>&1"


if WRITE_IGV_BED_ASSEMBLY:

    rule chimera_assembly_igv_bed:
        # BED track for IGV: one row per candidate, spanning the specific
        # TE-overlapping exon (not the whole transcript -- see
        # classify_chimera_assembly.py), colored by chimera_type. Uses the
        # cross-referenced table when available so junction-confirmed
        # candidates get a higher score (easy to filter/sort on in IGV).
        input:
            candidates=_chimera_assembly_summary_input(),
        output:
            "results/chimera/transcript_evidence/igv/transcripts.bed",
        threads: get_resources("chimera_assembly_igv_bed")["threads"]
        resources:
            mem_mb=get_resources("chimera_assembly_igv_bed")["mem_mb"],
            runtime=get_resources("chimera_assembly_igv_bed")["runtime"],
        log:
            "results/pipeline_info/logs/chimera_assembly/igv_bed.log",
        shell:
            "python3 {SCRIPTS_DIR}/chimera_assembly_to_igv_bed.py "
            "--candidates {input.candidates} --out {output} > {log} 2>&1"


# -----------------------------------------------------------------------------
# Assembly sample-QC: the PCA / sample-distance view the read screen already
# has, over this screen's own per-sample data (tpm_matrix.tsv.gz).
#
# These live HERE rather than in sample_qc.smk on purpose: that file is
# included only when the JUNCTION screen is on (Snakefile), and guard 27 pins
# that the assembly screen runs independently of it. Putting them there would
# silently drop this view whenever assembly runs alone.
#
# transform is log2, not vst/rlog: the matrix is already TPM, so DESeq2's
# count-based normalization does not apply. Thresholds are borrowed from
# chimera.junction.qc rather than adding a parallel config block -- the two
# views answer the same question and there is no evidence they want different
# cut-offs. Split them if that ever stops being true.
# -----------------------------------------------------------------------------
rule chimera_assembly_qc_transform:
    input:
        tpm="results/chimera/transcript_evidence/tpm_matrix.tsv.gz",
    output:
        "results/chimera/qc/assembly_log2_counts.tsv.gz",
    params:
        samples=config["samples"],
        min_samples_present=CHIMERA_QC["min_samples_present"],
        min_total_counts=CHIMERA_QC["min_total_counts"],
    threads: get_resources("chimera_assembly_qc_transform")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_assembly_qc_transform"),
        runtime=get_resources("chimera_assembly_qc_transform")["runtime"],
    log:
        "results/pipeline_info/logs/chimera_assembly/qc_transform.log",
    conda:
        CHIMERA_QC_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/sample_qc.R "
        "--transform assembly {input.tpm} {params.samples} log2 "
        "{params.min_samples_present} {params.min_total_counts} "
        "{output} > {log} 2>&1"


rule chimera_assembly_qc:
    input:
        transformed="results/chimera/qc/assembly_log2_counts.tsv.gz",
    output:
        pca="results/chimera/qc/assembly_pca_log2_mqc.json",
        heatmap="results/chimera/qc/assembly_heatmap_log2_mqc.json",
    params:
        samples=config["samples"],
        min_events=CHIMERA_QC["min_events"],
    threads: get_resources("chimera_assembly_qc")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("chimera_assembly_qc"),
        runtime=get_resources("chimera_assembly_qc")["runtime"],
    log:
        "results/pipeline_info/logs/chimera_assembly/qc_plots.log",
    conda:
        CHIMERA_QC_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/sample_qc.R "
        "--plots assembly {input.transformed} {params.samples} "
        "{params.min_events} log2 {output.pca} {output.heatmap} > {log} 2>&1"
