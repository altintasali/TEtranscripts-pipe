# Contrasts, the all_*() target lists, benchmark bookkeeping.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

def all_tecount_tables():
    return expand("results/tecount/{sample}.cntTable.gz", sample=SAMPLES)


def all_telocal_outputs():
    """TElocal per-sample count tables, counts matrix + sample-QC artifacts,
    and summary barplots for the `all` target (Snakefile). The QC view only
    runs when telocal.qc.enabled is true; the summary barplots always
    render."""
    files = expand("results/telocal/{sample}.cntTable.gz", sample=SAMPLES)
    if TELOCAL_QC_ENABLED:
        transform = TELOCAL_QC["pca_transform"]
        files += [
            "results/telocal/counts_matrix.tsv.gz",
            "results/telocal/qc/counts_matrix.tsv.gz",
            f"results/telocal/qc/{transform}_counts.tsv.gz",
            f"results/telocal/qc/pca_{transform}_mqc.json",
            f"results/telocal/qc/heatmap_{transform}_mqc.json",
        ]
    files += [
        "results/telocal/qc/telocal_assignment_mqc.json",
        "results/telocal/qc/telocal_te_class_mqc.json",
        "results/telocal/telocal_locations.bed",
    ]
    return files


def all_trim_outputs():
    """One trimmed fastq path per sample (only when trimming is enabled),
    used by the `trimming_only` convenience target (Snakefile)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1.fq.gz"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed.fq.gz"
        )
        for s in SAMPLES
    ]


def all_fastqc_reports():
    """One FastQC .zip report path per trimmed sample (only when trimming is
    enabled), used by the multiqc rule so it scans results/trimming/ for the
    FastQC + TrimGalore! reports without depending on the trimmed fastq files
    themselves (which may be temp()-deleted after alignment)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1_fastqc.zip"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed_fastqc.zip"
        )
        for s in SAMPLES
    ]


def all_raw_fastqc_reports():
    """FastQC .zip report paths for the raw (merged) input fastqs, one per
    sample/read. Run unconditionally so the MultiQC report always covers the
    untrimmed input, regardless of the optional `trimming` step."""
    return [
        f"results/fastqc/raw/{s}_R{read}_fastqc.zip"
        for s in SAMPLES
        for read in (1, 2)
        if read == 1 or _is_paired(s)
    ]


def all_benchmark_files():
    """Every benchmark file this configuration will produce, so the
    benchmark_summary rule (qc.smk) aggregates exactly the rules that ran
    into the MultiQC resource-usage section. Kept in lockstep with the
    `benchmark:` declarations in the rules -- only the conditional ones
    (merging, trimming, gunzip, RSeQC bed12 conversion, strandedness
    auto-detection) need to be gated here. The multiqc and
    benchmark_summary rules' own benchmarks are deliberately excluded to
    avoid a cyclic dependency (their resource use is negligible)."""
    B = "results/pipeline_info/benchmarks"
    files = [
        "results/pipeline_info/benchmarks/software_versions/software_versions.txt",
        "results/pipeline_info/benchmarks/gene_name_lookup/gene_name_lookup.txt",
    ]
    # The star_index rule only exists when we build the index ourselves; with
    # star.build_index: false the external index is a plain input and no such
    # job (or benchmark) is ever produced.
    if STAR_BUILD_INDEX:
        files.append(
            "results/pipeline_info/benchmarks/star_index/star_index.txt"
        )
    if STRAND_CHECK_SAMPLES:
        files.append(
            "results/pipeline_info/benchmarks/strandedness_check/"
            "strandedness_check.txt"
        )
    # The BED12 gene model is built on every run: rseqc_read_distribution and
    # rseqc_gene_body_coverage (bam_qc.smk) consume it for ALL samples, not
    # just the strandedness ones. It used to be gated on AUTO_SAMPLES, from
    # when RSeQC strandedness was its only consumer.
    files += [
        "results/pipeline_info/benchmarks/gtf_to_genepred/gtf_to_genepred.txt",
        "results/pipeline_info/benchmarks/genepred_to_bed12/genepred_to_bed12.txt",
    ]
    for stem in REFERENCE_GZ_SOURCES:
        files.append(f"results/pipeline_info/benchmarks/gunzip_reference/{stem}.txt")
    for s in SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/star_align/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_sort/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_index/{s}.txt",
            f"results/pipeline_info/benchmarks/tecount/{s}.txt",
            f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R1.txt",
            # bam_qc.smk -- always on, for every sample. Omitted here until
            # now, so the table silently hid them; gene-body coverage in
            # particular is one of the slowest per-sample steps.
            f"results/pipeline_info/benchmarks/samtools_flagstat/{s}.txt",
            f"results/pipeline_info/benchmarks/rseqc_read_distribution/{s}.txt",
            f"results/pipeline_info/benchmarks/rseqc_gene_body_coverage/{s}.txt",
        ]
        if _is_paired(s):
            files.append(f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R2.txt")
        # Lane-merging (cat_fastq) only for samples with multiple lanes/runs.
        for read in (1, 2):
            if (read == 1 or _is_paired(s)) and len(sample_fastqs(s, read)) > 1:
                files.append(
                    f"results/pipeline_info/benchmarks/cat_fastq/{s}_R{read}.txt"
                )
        if TRIM_ENABLED:
            files.append(
                f"results/pipeline_info/benchmarks/"
                f"{'trim_galore_pe' if _is_paired(s) else 'trim_galore_se'}/{s}.txt"
            )
    for s in STRAND_CHECK_SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/rseqc_infer_experiment/{s}.txt",
            f"results/pipeline_info/benchmarks/determine_strandedness/{s}.txt",
        ]
    # Chimera-screen rules only run when the chimera stage is enabled.
    if CHIMERA_CHIMERIC_READS_ENABLED:
        files += [
            "results/pipeline_info/benchmarks/annotation_to_bed/annotation_to_bed.txt",
            "results/pipeline_info/benchmarks/chimera_chimeric_reads_counts/chimera_chimeric_reads_counts.txt",
            "results/pipeline_info/benchmarks/chimera_chimeric_reads_highlights/"
            "chimera_chimeric_reads_highlights.txt",
            "results/pipeline_info/benchmarks/chimera_evidence/"
            "chimera_evidence.txt",
            "results/pipeline_info/benchmarks/chimera_evidence_guide/"
            "chimera_evidence_guide.txt",
            "results/pipeline_info/benchmarks/chimera_candidates_table/"
            "chimera_candidates_table.txt",
            "results/pipeline_info/benchmarks/chimera_chimeric_reads_te_type/"
            "chimera_chimeric_reads_te_type.txt",
        ]
        files.append(
            f"{B}/chimera_chimeric_reads_qc_barplot/chimera_chimeric_reads_qc_barplot.txt"
        )
        if TELOCAL_ENABLED:
            # The reads screen cross-references TElocal only when it ran.
            files.append(f"{B}/chimera_telocal_index/chimera_telocal_index.txt")
        for s in SAMPLES:
            files += [
                f"{B}/star_filter_primary/{s}.txt",
                f"{B}/chimera_chimeric_reads_classify/{s}.txt",
                f"{B}/chimera_chimeric_reads_qc/{s}.txt",
            ]
            if TELOCAL_ENABLED:
                files.append(f"{B}/chimera_telocal_annotate/{s}.txt")
            if config["chimera"]["chimeric_reads"]["outputs"]["write_igv_bed"]:
                files.append(
                    f"results/pipeline_info/benchmarks/chimera_chimeric_reads_igv_bed/{s}.txt"
                )
        if (config["chimera"]["chimeric_reads"]["outputs"]["write_counts_matrix"]
                and config["chimera"]["chimeric_reads"]["qc"].get("enabled", False)):
            transform = config["chimera"]["chimeric_reads"]["qc"]["pca_transform"]
            files += [
                f"results/pipeline_info/benchmarks/"
                f"chimera_chimeric_reads_sample_qc_transform/{transform}.txt",
                f"results/pipeline_info/benchmarks/chimera_chimeric_reads_sample_qc/{transform}.txt",
            ]
    # Assembly screen: its own STAR pass, StringTie, and everything after.
    # None of this was listed, so the second-heaviest stage in the workflow
    # was invisible in the resource table.
    if CHIMERA_ASSEMBLY_ENABLED:
        files += [
            f"{B}/stringtie_merge/stringtie_merge.txt",
            f"{B}/chimera_assembly_classify/chimera_assembly_classify.txt",
            f"{B}/chimera_assembly_quantify/chimera_assembly_quantify.txt",
            f"{B}/chimera_assembly_summary_mqc/chimera_assembly_summary_mqc.txt",
        ]
        for s in SAMPLES:
            files += [
                f"{B}/star_align_for_assembly/{s}.txt",
                f"{B}/stringtie_assemble/{s}.txt",
                f"{B}/stringtie_requantify/{s}.txt",
            ]
            if KEEP_ASSEMBLY_BAM:
                files.append(f"{B}/chimera_assembly_bam_index/{s}.txt")
        if CHIMERA_CHIMERIC_READS_ENABLED:
            files.append(
                f"{B}/chimera_assembly_cross_evidence/chimera_assembly_cross_evidence.txt"
            )
        if config["chimera"]["assembly"]["outputs"]["write_igv_bed"]:
            files.append(f"{B}/chimera_assembly_igv_bed/chimera_assembly_igv_bed.txt")
        if config["chimera"]["assembly"]["outputs"]["write_gene_te_chimera_counts"]:
            files.append(
                f"{B}/chimera_assembly_aggregate_counts/chimera_assembly_aggregate_counts.txt"
            )
        # The assembly QC view is log2-fixed (its outputs carry no
        # {transform} wildcard), unlike the reads screen's. Borrows
        # chimera.chimeric_reads.qc.enabled (no qc block of its own).
        if config["chimera"]["chimeric_reads"]["qc"].get("enabled", False):
            files += [
                f"{B}/chimera_assembly_qc_transform/chimera_assembly_qc_transform.txt",
                f"{B}/chimera_assembly_qc/chimera_assembly_qc.txt",
            ]

    # Cohort-wide STAR 2-pass: a pass-1 alignment per sample plus one merge.
    if STAR_TWO_PASS == "cohort":
        files.append(f"{B}/star_merge_junctions/merge.txt")
        for s in SAMPLES:
            files.append(f"{B}/star_align_pass1/{s}.txt")

    # SJ.out.tab junction screen: no extra STAR pass (reuses star_align's
    # own SJ.out.tab), so only its own classify/counts/QC rules are new.
    if CHIMERA_SPLICE_JUNCTIONS_ENABLED:
        files.append(f"{B}/chimera_splice_junctions_counts/chimera_splice_junctions_counts.txt")
        files.append(
            f"{B}/chimera_splice_junctions_summary_mqc/"
            "chimera_splice_junctions_summary_mqc.txt"
        )
        for s in SAMPLES:
            files.append(f"{B}/chimera_splice_junctions_classify/{s}.txt")
            if config["chimera"]["splice_junctions"]["outputs"]["write_igv_bed"]:
                files.append(
                    f"results/pipeline_info/benchmarks/chimera_splice_junctions_igv_bed/{s}.txt"
                )
        if (config["chimera"]["splice_junctions"]["outputs"]["write_counts_matrix"]
                and config["chimera"]["splice_junctions"]["outputs"]["write_gene_te_chimera_counts"]):
            files.append(
                f"{B}/chimera_splice_junctions_aggregate_counts/"
                "chimera_splice_junctions_aggregate_counts.txt"
            )
        if (config["chimera"]["splice_junctions"]["outputs"]["write_counts_matrix"]
                and config["chimera"]["splice_junctions"]["qc"].get("enabled", False)):
            transform = config["chimera"]["splice_junctions"]["qc"]["pca_transform"]
            files += [
                f"{B}/chimera_splice_junctions_qc_transform/{transform}.txt",
                f"{B}/chimera_splice_junctions_qc/{transform}.txt",
            ]

    # Report-assembly rules that are siblings of benchmark_summary (they do
    # not depend on it, so listing them just orders them earlier).
    files += [
        f"{B}/config_used/config_used.txt",
        f"{B}/evidence_overview/evidence_overview.txt",
    ]

    # TEcounts sample-QC rules only run when tetranscripts.qc.enabled.
    if TECOUNT_QC_ENABLED:
        transform = TECOUNT_QC["pca_transform"]
        files += [
            "results/pipeline_info/benchmarks/tecount_counts/tecount_counts.txt",
            f"results/pipeline_info/benchmarks/tecount_qc_transform/{transform}.txt",
            f"results/pipeline_info/benchmarks/tecount_qc/{transform}.txt",
        ]
    # The summary-barplot rule runs on every run (raw cntTables only).
    files.append(
        "results/pipeline_info/benchmarks/tecount_summary/tecount_summary.txt"
    )
    # TElocal rules only run when the telocal stage is enabled; its QC view
    # additionally requires telocal.qc.enabled.
    if TELOCAL_ENABLED:
        files += [
            "results/pipeline_info/benchmarks/telocal_locind/locind.txt",
            "results/pipeline_info/benchmarks/telocal_summary/telocal_summary.txt",
            "results/pipeline_info/benchmarks/telocal_locations/locations.txt",
        ]
        for s in SAMPLES:
            files.append(f"results/pipeline_info/benchmarks/telocal/{s}.txt")
        if TELOCAL_QC_ENABLED:
            transform = TELOCAL_QC["pca_transform"]
            files += [
                "results/pipeline_info/benchmarks/telocal_counts/telocal_counts.txt",
                f"results/pipeline_info/benchmarks/telocal_qc_transform/{transform}.txt",
                f"results/pipeline_info/benchmarks/telocal_qc/{transform}.txt",
            ]
    return sorted(set(files))


def allocated_resources_by_rule():
    """{rule: {"threads", "mem_mb"}} for every rule that has benchmark files
    (the benchmark_summary rule's input), read from resources.yaml -- the
    per-job allocation against which the benchmark_summary script computes
    CPU/RAM efficiency. Iterated in Snakemake's own rule order (workflow.rules,
    an OrderedDict) so the resource-usage table lists rules in the order they
    appear in the workflow, not alphabetically."""
    benchmark_rules = {
        os.path.basename(os.path.dirname(path)) for path in all_benchmark_files()
    }
    out = {}
    for r in workflow.rules:
        if r.name in benchmark_rules:
            out[r.name] = get_resources(r.name)
    return out
