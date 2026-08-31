# Generated per-tool conda environments.
#
# Part of the former single common.smk (1,217 lines doing eight jobs).
# Included by rules/common.smk in a fixed order -- these files are NOT
# independent: each builds on names the previous ones defined.

# -----------------------------------------------------------------------------
# Conda environments for every tool in the workflow, generated from
# config["versions"]. Editing a version string there and re-running is all
# that's needed -- Snakemake/conda resolves and downloads the matching build
# automatically (via `--sdm conda` / `--use-conda`), no manual install
# required. This intentionally does NOT rely on snakemake-wrapper-bundled
# environments, because a wrapper's pinned tag moves all of its tools'
# versions together -- it can't give you e.g. STAR 2.7.11b + samtools 1.21
# independently if that particular wrapper tag happens to pin a different
# samtools. Generating one small env per tool from explicit version strings
# is the only way to match an arbitrary external reference (e.g. a specific
# nf-core/rnaseq release's tool versions) exactly.
# -----------------------------------------------------------------------------
# Absolute path: rules in this repo are included from workflow/Snakefile, and
# Snakemake resolves relative paths declared inside an included file against
# that file's directory (workflow/rules/), not the run directory -- so a
# relative "workflow/envs/generated" here would make rule env: references point
# at workflow/rules/workflow/envs/generated/ (which never exists). The env
# files are written at parse time relative to the process CWD, so an absolute
# path is the one form both sides agree on.
GENERATED_ENV_DIR = os.path.abspath("workflow/envs/generated")
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)


def _write_env(name, dependencies, pip_dependencies=None):
    dependencies = list(dependencies)
    if pip_dependencies:
        dependencies.append({"pip": list(pip_dependencies)})
    path = f"{GENERATED_ENV_DIR}/{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(
            {"channels": ["bioconda", "conda-forge"], "dependencies": dependencies},
            fh,
            sort_keys=False,
        )
    return path


STAR_ENV = _write_env("star", [f"star={V['star']}"])
SAMTOOLS_ENV = _write_env("samtools", [f"samtools={V['samtools']}"])
STRINGTIE_ENV = _write_env("stringtie", [f"stringtie={V['stringtie']}"])
# trim-galore brings cutadapt (its core dependency) along automatically;
# fastqc is added explicitly because trim_galore's --fastqc_args (nf-core/
# rnaseq default) shells out to it, and its reports feed the MultiQC report.
TRIM_GALORE_ENV = _write_env(
    "trim_galore", [f"trim-galore={V['trim_galore']}", f"fastqc={V['fastqc']}"]
)
# Standalone FastQC env for the always-on raw FastQC (qc.smk), independent of
# the optional trimming step.
FASTQC_ENV = _write_env("fastqc", [f"fastqc={V['fastqc']}"])
# python>=3.9 floor: without it, conda's solver can backtrack all the way to
# an ancient RSeQC/MultiQC build (seen in practice: Python 3.6, from ~2021)
# to find something that resolves at all -- and those old builds pull in a
# pysam/htslib linked against OpenSSL 1.0, which doesn't exist on modern
# systems (`ImportError: libcrypto.so.1.0.0: cannot open shared object
# file`). Pinning a modern floor forces the solver toward current,
# self-consistent builds instead.
RSEQC_ENV = _write_env("rseqc", [f"rseqc={V['rseqc']}", "python>=3.9"])
MULTIQC_ENV = _write_env("multiqc", [f"multiqc={V['multiqc']}", "python>=3.9"])

# Absolute path to workflow/scripts: shell directives that run the workflow's
# own python/R scripts need a path that resolves the same way from the run
# directory (the repo root) as the included rules file is parsed from.
SCRIPTS_DIR = os.path.abspath("workflow/scripts")

# The MultiQC custom config, resolved the same way. Declared as a rule INPUT
# (not just interpolated into the shell command) so editing it re-runs the
# report: it controls section order, module naming and version detection, so
# a change to it changes the output, and an untracked change silently did
# not. Excluded from the rule's search directories -- see the multiqc rule.
MULTIQC_CONFIG = os.path.abspath("workflow/default-config/multiqc_config.yaml")

# Chimera sample-QC (PCA / sample clustering) runs in R with DESeq2
# (nf-core/rnaseq style); deseq2 + r-base come from conda. A python>=3.9
# floor guards the solver against ancient, mutually-incompatible builds.
# Only referenced when the chimera stage is enabled.
CHIMERA_QC_ENV = _write_env(
    "chimera_qc",
    [
        f"bioconductor-deseq2={V['deseq2']}",
        f"r-base={V['r_base']}",
        "r-pheatmap",
        "python>=3.9",
    ],
)

# TEtranscripts is installed from PyPI rather than bioconda: the bioconda
# recipe's run dependencies pin an ancient bioconductor-deseq (DESeq v1),
# which is not used at runtime (only DESeq2 is) and can only coexist with
# R 4.0-era packages -- so the conda package is unsolvable together with a
# modern bioconductor-deseq2/r-base on any platform. The PyPI package is
# pure Python (depends only on pysam); DESeq2 and R are provided by conda
# as before, so the deseq2/r_base version pins still apply.
TETRANSCRIPTS_ENV = _write_env(
    "tetranscripts",
    [
        f"bioconductor-deseq2={V['deseq2']}",
        f"r-base={V['r_base']}",
        "pysam",
        "pip",
    ],
    pip_dependencies=[f"TEtranscripts=={V['tetranscripts']}"],
)

UCSC_TOOLS_ENV = _write_env(
    "ucsc_tools",
    [
        f"ucsc-gtftogenepred={V['ucsc_gtftogenepred']}",
        f"ucsc-genepredtobed={V['ucsc_genepredtobed']}",
    ],
)
