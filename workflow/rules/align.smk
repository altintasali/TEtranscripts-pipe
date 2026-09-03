import os


def _chim_star_args(wildcards):
    """STAR chimeric-alignment parameters from config["chimera"]["reads"]["star"], or
    nothing when the chimera screen is disabled (chimera.enabled: false).
    --chimOutType stays fixed (Junctions + WithinBAM SoftClip): the junction
    file feeds the chimera screen and the SA tags embedded in the BAM let the
    same alignment be re-inspected in IGV."""
    if not CHIMERA_READS_ENABLED:
        return ""
    c = config["chimera"]["reads"]["star"]
    return " ".join(
        [
            f"--chimSegmentMin {c['segment_min']}",
            f"--chimJunctionOverhangMin {c['overhang_min']}",
            f"--chimScoreDropMax {c['score_drop_max']}",
            str(c.get("extra", "")),
        ]
    ).strip()


rule star_align:
    # Aligns one sample against the shared genome index. Output BAM is left
    # "Unsorted" (STAR's native read order) on purpose: TEtranscripts/TEcount
    # require either unsorted or queryname-sorted input, see
    # https://github.com/mhammell-laboratory/TEtranscripts#recommendations-for-tetranscripts-input-files
    # Runs the STAR version pinned in config["versions"]["star"].
    # When the chimera stage is enabled (config["chimera"]["reads"]["enabled"]), the
    # same run also detects chimeric junctions (--chimOutType Junctions
    # WithinBAM SoftClip) -- the {sample}_Chimeric.out.junction file that the
    # chimera screen consumes is a side output of this same alignment, no
    # separate alignment step.
    #
    # config star.two_pass controls STAR 2-pass mapping (see common.smk and
    # the README): "none" (default) adds nothing; "per_sample" is STAR's own
    # --twopassMode Basic; "cohort" additionally depends on
    # results/star_pass1/merged_SJ.out.tab (rules/star_two_pass.smk), so this
    # rule doesn't start for ANY sample until every sample's pass-1 has run.
    #
    # STAR's scratch directory (default: prefix/_STARtmp) is redirected to
    # --outTmpDir under config["star"]["tmpdir"] (default: the OS temp dir)
    # so it never clutters results/ -- it holds ~a BAM's worth of per-thread
    # chunks that STAR merges and then deletes. A trap removes it even when
    # STAR crashes before its own cleanup (the known benign exit-time crash,
    # see below), so nothing lingers in /tmp either.
    input:
        unpack(star_input),
        idx=STAR_INDEX_DIR,
        # merged_SJ.out.tab (cohort-wide 2-pass only) is what makes this
        # rule wait on every sample's star_align_pass1 + the
        # star_merge_junctions aggregation -- see rules/star_two_pass.smk.
        **({"merged_sj": "results/star_pass1/merged_SJ.out.tab"}
           if STAR_TWO_PASS == "cohort" else {}),
    output:
        aln="results/star/{sample}_Aligned.out.bam",
        log_final="results/star/{sample}_Log.final.out",
        sj="results/star/{sample}_SJ.out.tab",
        # The {sample}_Chimeric.out.junction file (the chimera screen's
        # input) is only declared when the chimera stage is enabled;
        # otherwise it is not produced and not required.
        **({"chim": "results/star/{sample}_Chimeric.out.junction"}
           if CHIMERA_READS_ENABLED else {}),
    params:
        reads=star_reads_param,
        read_command=star_read_command_param,
        prefix=lambda wc: f"results/star/{wc.sample}_",
        tmpdir=lambda wc: os.path.join(STAR_TMPDIR, f"star_{wc.sample}"),
        chim=_chim_star_args,
        chim_out_type=(
            "--chimOutType Junctions WithinBAM SoftClip " if CHIMERA_READS_ENABLED else ""
        ),
        # STAR 2-pass mapping (config star.two_pass, see common.smk and the
        # README): "per_sample" is STAR's own --twopassMode Basic;
        # "cohort" points STAR at the pooled novel-junction file built by
        # star_merge_junctions (STAR accepts SJ.out.tab-formatted files
        # directly via --sjdbFileChrStartEnd, no conversion needed -- STAR
        # manual section 9.1); "none" adds nothing.
        two_pass=lambda wc, input: (
            "--twopassMode Basic"
            if STAR_TWO_PASS == "per_sample"
            else f"--sjdbFileChrStartEnd {input.merged_sj}"
            if STAR_TWO_PASS == "cohort"
            else ""
        ),
        extra=config["star"]["extra"],
    threads: get_resources("star_align")["threads"]
    resources:
        mem_mb=get_resources("star_align")["mem_mb"],
        runtime=get_resources("star_align")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/star_align/{sample}.txt",
    log:
        "results/pipeline_info/logs/star/align/{sample}.log",
    conda:
        STAR_ENV
    shell:
        # See the star_index rule's comment: STAR has a known benign
        # crash-on-exit bug (confirmed by its author) that can happen
        # after alignment has already finished and all output files were
        # written -- so a non-zero exit is checked against actual output
        # completeness rather than trusted outright. --outTmpDir must be
        # ABSENT before STAR starts: STAR creates it itself and aborts with
        # "could not make temporary directory" if it already exists, so it is
        # rm -rf'd beforehand to clear stale leftovers from a crashed job.
        "mkdir -p results/star && "
        "rm -rf {params.tmpdir} && "
        "trap 'rm -rf {params.tmpdir}' EXIT; "
        "(STAR --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {params.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outTmpDir {params.tmpdir} "
        "--outSAMtype BAM Unsorted "
        "{params.chim_out_type}"
        "{params.chim} "
        "{params.two_pass} "
        "{params.extra} "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether alignment output "
        "is actually complete anyway (see rule comment re: known benign "
        "STAR exit-time crash)' >> {log}; "
        "if test -s {output.aln} && test -s {output.log_final} && "
        "grep -q 'ALL DONE!' {output.log_final}; then :; else "
        "echo '-- STAR failed for real; log tail --' >&2; "
        "tail -n 60 {log} >&2; exit 1; fi))"


def _samtools_sort_mem(wildcards):
    """samtools sort -m flag (max memory *per thread*), derived from this
    rule's mem_mb/threads resources so it scales with input/resources.yaml
    instead of being hardcoded. Total = -m x threads, so per-thread is capped
    at a fraction of the per-thread share to leave headroom for the sort's
    I/O buffers and keep peak usage inside the job's mem_mb budget.
    """
    res = get_resources("samtools_sort")
    threads = max(res["threads"], 1)
    per_thread = max(int(res["mem_mb"] * 0.8 // threads), 64)
    return f"-m {per_thread}M"


rule samtools_sort:
    # A coordinate-sorted + indexed copy is only needed for RSeQC/QC purposes;
    # TEcount/TEtranscripts always consume the unsorted STAR output above.
    # Runs the samtools version pinned in config["versions"]["samtools"].
    input:
        "results/star/{sample}_Aligned.out.bam",
    output:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam",
    params:
        extra=_samtools_sort_mem,
    threads: get_resources("samtools_sort")["threads"]
    resources:
        mem_mb=get_resources("samtools_sort")["mem_mb"],
        runtime=get_resources("samtools_sort")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/samtools_sort/{sample}.txt",
    log:
        "results/pipeline_info/logs/samtools/sort/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools sort {params.extra} -@ {threads} "
        "-o {output} {input} > {log} 2>&1"


rule samtools_index:
    input:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    params:
        extra="",
    threads: get_resources("samtools_index")["threads"]
    resources:
        mem_mb=get_resources("samtools_index")["mem_mb"],
        runtime=get_resources("samtools_index")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/samtools_index/{sample}.txt",
    log:
        "results/pipeline_info/logs/samtools/index/{sample}.log",
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools index {params.extra} -@ {threads} {input} {output} > {log} 2>&1"
