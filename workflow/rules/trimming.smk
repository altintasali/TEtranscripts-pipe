# -----------------------------------------------------------------------------
# Optional TrimGalore! adapter/quality trimming before STAR, mirroring the
# nf-core/rnaseq defaults (Trim Galore! is nf-core's default trimmer, run
# with --paired/--gzip and FastQC executed inside via --fastqc_args).
# Controlled by the `trimming` section of config.yaml; set enabled: false to
# skip this step entirely (STAR then reads the merged/raw fastqs directly).
#
# Two rules because Snakemake requires a rule's output patterns to be static
# and TrimGalore! names paired vs single-end outputs differently. Only one of
# them ever triggers per sample (paired vs single-end is a per-sample
# property), so they can share the shell and helper functions below.
# -----------------------------------------------------------------------------


def _trim_inputs(wildcards):
    d = {"r1": merged_fastq_path(wildcards.sample, 1)}
    if _is_paired(wildcards.sample):
        d["r2"] = merged_fastq_path(wildcards.sample, 2)
    return d


def _trim_reads_param(wildcards):
    d = _trim_inputs(wildcards)
    return d["r1"] + (" " + d["r2"] if "r2" in d else "")


def _trim_paired_flag(wildcards):
    return "--paired" if _is_paired(wildcards.sample) else ""


def _nextseq_param(wildcards):
    n = int(config["trimming"].get("trim_nextseq", 0) or 0)
    return f"--nextseq {n}" if n > 0 else ""


def _trim_outdir(wildcards):
    return "results/trimming"


_TRIM_SHELL = (
    # TrimGalore! records the INPUT FILE NAME inside its trimming report
    # ("Input filename: ..."), and MultiQC's cutadapt module takes the sample
    # name from that line -- not from the report's filename. --basename fixes
    # the trimmed fastqs and their FastQC zips but has no effect on it.
    #
    # A single-lane sample is handed its RAW fastq path (see
    # merged_fastq_path), so on real data the report recorded the sequencing
    # facility's name and MultiQC grew an extra General Statistics row per
    # sample -- carrying only "Trimmed bases", never merging with that
    # sample's other rows. Multi-lane samples were unaffected (their input is
    # already {sample}_R{read}.fastq.gz), which is why the test data never
    # showed it.
    #
    # So the input is symlinked to a canonical {sample}_R{n} name first and
    # TrimGalore! is run on the symlink -- the same fix fastqc_raw uses, and
    # for the same reason: renaming the output cannot correct a name the tool
    # has already written inside the file. The _R1/_R2 suffix is in
    # extra_fn_clean_exts, so MultiQC resolves it to {sample}.
    "mkdir -p {params.outdir} && "
    "stage={params.outdir}/.stage_{wildcards.sample} && "
    "rm -rf \"$stage\" && mkdir -p \"$stage\" && "
    "trap 'rm -rf \"$stage\"' EXIT && "
    # Positional parameters rather than a string accumulator: TrimGalore!
    # aborts on a filename list with leading whitespace, which is exactly what
    # "$acc $next" produces on the first iteration.
    "set --; i=1; "
    "for f in {params.reads}; do "
    # keep the original compression: a .gz symlinked as plain (or vice
    # versa) makes TrimGalore! misread the stream.
    "  case \"$f\" in *.gz) ext=.fastq.gz;; *) ext=.fastq;; esac; "
    "  ln -s \"$(cd \"$(dirname \"$f\")\" && pwd)/$(basename \"$f\")\" "
    "    \"$stage/{wildcards.sample}_R${{i}}$ext\"; "
    "  set -- \"$@\" \"$stage/{wildcards.sample}_R${{i}}$ext\"; "
    "  i=$((i+1)); "
    "done && "
    "trim_galore {params.paired} --gzip --cores {threads} "
    "{params.nextseq} {params.extra} "
    "--fastqc_args '-t {threads}' "
    "--basename {wildcards.sample} --output_dir {params.outdir} "
    "\"$@\" "
    "> {log} 2>&1 && "
    "python3 workflow/scripts/check_nonempty_fastq.py "
    "{params.outdir}/{wildcards.sample}_val_1.fq.gz "
    "{params.outdir}/{wildcards.sample}_val_2.fq.gz "
    "{params.outdir}/{wildcards.sample}_trimmed.fq.gz"
)


rule trim_galore_pe:
    # Paired-end: --basename fixes the output names ({sample}_val_1/2.fq.gz)
    # so they don't depend on the input file's basename -- the merged path
    # (results/fastq/) and a raw single-lane path differ, but both must map
    # to the same outputs. The FastQC reports (run inside TrimGalore! via
    # --fastqc_args) are declared as outputs too: the multiqc rule depends on
    # them rather than on the (possibly temp()-deleted) trimmed fastqs.
    input:
        unpack(_trim_inputs),
    output:
        fq1=_maybe_temp("results/trimming/{sample}_val_1.fq.gz", KEEP_TRIMMED_FASTQ),
        fq2=_maybe_temp("results/trimming/{sample}_val_2.fq.gz", KEEP_TRIMMED_FASTQ),
        fq1_fastqc="results/trimming/{sample}_val_1_fastqc.zip",
        fq2_fastqc="results/trimming/{sample}_val_2_fastqc.zip",
    params:
        paired=_trim_paired_flag,
        reads=_trim_reads_param,
        nextseq=_nextseq_param,
        extra=config["trimming"].get("extra", ""),
        outdir=_trim_outdir,
    threads: get_resources("trim_galore_pe")["threads"]
    resources:
        mem_mb=get_resources("trim_galore_pe")["mem_mb"],
        runtime=get_resources("trim_galore_pe")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/trim_galore_pe/{sample}.txt",
    log:
        "results/pipeline_info/logs/trimming/{sample}.log",
    conda:
        TRIM_GALORE_ENV
    shell:
        _TRIM_SHELL


rule trim_galore_se:
    # Single-end: TrimGalore! names the output {sample}_trimmed.fq.gz.
    input:
        unpack(_trim_inputs),
    output:
        fq=_maybe_temp("results/trimming/{sample}_trimmed.fq.gz", KEEP_TRIMMED_FASTQ),
        fastqc="results/trimming/{sample}_trimmed_fastqc.zip",
    params:
        paired=_trim_paired_flag,
        reads=_trim_reads_param,
        nextseq=_nextseq_param,
        extra=config["trimming"].get("extra", ""),
        outdir=_trim_outdir,
    threads: get_resources("trim_galore_se")["threads"]
    resources:
        mem_mb=get_resources("trim_galore_se")["mem_mb"],
        runtime=get_resources("trim_galore_se")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/trim_galore_se/{sample}.txt",
    log:
        "results/pipeline_info/logs/trimming/{sample}.se.log",
    conda:
        TRIM_GALORE_ENV
    shell:
        _TRIM_SHELL
