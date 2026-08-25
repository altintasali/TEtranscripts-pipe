import os


rule star_align_pass1:
    # Cohort-wide 2-pass mapping, stage 1 (STAR manual section 9.1): a
    # throwaway alignment per sample whose only purpose is discovering
    # splice junctions from this sample's own reads. --outSAMtype None skips
    # writing a BAM entirely (not used downstream) to save time/disk; the
    # genome index is the same annotation-guided STAR_INDEX_DIR used by the
    # real alignment, so ref.gtf's junctions are already covered without
    # extra flags here. No chimeric detection -- that only runs in the final
    # alignment (rules/align.smk's star_align).
    input:
        unpack(star_input),
        idx=STAR_INDEX_DIR,
    output:
        log_final="results/star_pass1/{sample}_Log.final.out",
        sj="results/star_pass1/{sample}_SJ.out.tab",
    params:
        reads=star_reads_param,
        read_command=star_read_command_param,
        prefix=lambda wc: f"results/star_pass1/{wc.sample}_",
        tmpdir=lambda wc: os.path.join(STAR_TMPDIR, f"star_pass1_{wc.sample}"),
    threads: get_resources("star_align_pass1")["threads"]
    resources:
        mem_mb=get_resources("star_align_pass1")["mem_mb"],
        runtime=get_resources("star_align_pass1")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/star_align_pass1/{sample}.txt",
    log:
        "results/pipeline_info/logs/star/align_pass1/{sample}.log",
    conda:
        STAR_ENV
    shell:
        # Same benign-exit-crash tolerance as star_align (rules/align.smk) --
        # see that rule's comment. Success here means SJ.out.tab exists and
        # Log.final.out reports completion; there is no BAM to check.
        "mkdir -p results/star_pass1 && "
        "rm -rf {params.tmpdir} && "
        "trap 'rm -rf {params.tmpdir}' EXIT; "
        "(STAR --runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--readFilesIn {params.reads} "
        "{params.read_command} "
        "--outFileNamePrefix {params.prefix} "
        "--outTmpDir {params.tmpdir} "
        "--outSAMtype None "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether pass-1 output is "
        "actually complete anyway (see star_align rule comment re: known "
        "benign STAR exit-time crash)' >> {log}; "
        "if test -e {output.sj} && test -s {output.log_final} && "
        "grep -q 'ALL DONE!' {output.log_final}; then :; else "
        "echo '-- STAR pass-1 failed for real; log tail --' >&2; "
        "tail -n 60 {log} >&2; exit 1; fi))"


rule star_merge_junctions:
    # Cohort-wide 2-pass mapping, stage 2 (STAR manual section 9.1 step 2):
    # pools every sample's pass-1 novel junctions into one filtered file,
    # fed to every sample's final alignment via --sjdbFileChrStartEnd
    # (rules/align.smk's star_align). This is the DAG sync point -- no
    # sample's final alignment can start until every sample's pass-1 has
    # finished. A junction supported by only a couple of reads in any single
    # sample, but recurrently so across many samples, survives here even
    # though no single sample's pass-1 alone would have kept it.
    input:
        sj=expand("results/star_pass1/{sample}_SJ.out.tab", sample=SAMPLES),
    output:
        merged="results/star_pass1/merged_SJ.out.tab",
    params:
        min_reads=config["star"].get("two_pass_min_unique_reads", 3),
    threads: get_resources("star_merge_junctions")["threads"]
    resources:
        mem_mb=get_scaled_mem_mb("star_merge_junctions"),
        runtime=get_resources("star_merge_junctions")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/star_merge_junctions/merge.txt",
    log:
        "results/pipeline_info/logs/star/merge_junctions.log",
    shell:
        "python3 {SCRIPTS_DIR}/merge_splice_junctions.py "
        "--sj-tables {input.sj} "
        "--min-unique-reads {params.min_reads} "
        "--out {output.merged} > {log} 2>&1"
