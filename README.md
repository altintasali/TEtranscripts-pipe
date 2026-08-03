# rnaseq-star-tetranscripts

A Snakemake workflow that quantifies genes **and** transposable elements (TEs)
from RNA-seq data with [TEtranscripts/TEcount](https://github.com/mhammell-laboratory/TEtranscripts),
aligning with STAR and auto-detecting library strandedness via RSeQC.

```
fastq --> STAR align --> [samtools sort/index --> RSeQC infer_experiment --> strandedness] --> TEcount (per sample)
                                                                                             --> TEtranscripts (per contrast, only if a "condition" column is present)
                                                                                             --> MultiQC
```

## Quick start

Requires conda (e.g. [Miniforge](https://github.com/conda-forge/miniforge) or
Miniconda; conda ≥23.10 already uses the fast libmamba solver, so mamba is
optional).

```bash
git clone https://github.com/altintasali/rnaseq-star-tetranscripts.git
cd rnaseq-star-tetranscripts
conda env create -f environment.yaml        # snakemake + every analysis tool, installed once
conda activate rnaseq-star-tetranscripts

# Smoke test on the bundled tiny synthetic dataset (no real genome/reads needed):
snakemake --configfile config/test.yaml --cores 4

# Run on your own data (edit config/config.yaml + config/samples.csv first):
snakemake --cores 16
```

All tools (STAR, samtools, RSeQC, MultiQC, TEtranscripts, DESeq2, UCSC tools)
live directly in this environment, so no `--sdm conda` is needed — nothing to
download or solve per run.

## Configuration

You only need to touch two files to run an analysis:

**`config/config.yaml`** — reference paths and tool options:

| key | meaning |
|-----|---------|
| `ref.fasta` | genome FASTA (Ensembl/GENCODE). Plain or gzipped. |
| `ref.gtf` | gene annotation GTF. Plain or gzipped. |
| `ref.te_gtf` | **curated** TE GTF from the TEtranscripts authors (download the file matching your genome build from [mghlab.org/software/tetranscripts](https://www.mghlab.org/software/tetranscripts) — a generic RepeatMasker GTF will *not* work). Plain or gzipped. |
| `ref.sjdb_overhang` | `auto` (default) detects `max(read length) - 1` from your fastq files; set an integer to pin it. |
| `star.index` | where to (re)build the STAR index. |
| `star.extra` | alignment flags; pre-set to the TEtranscripts authors' multi-mapper recommendations. |
| `strandedness.min_fraction` | confidence threshold for RSeQC auto-detection. |
| `tetranscripts.*` | TEcount/TEtranscripts options (mode, padj, foldchange...). |

**`config/samples.csv`** — the design file. Columns in this exact order
(an nf-core/rnaseq samplesheet works as-is):

| column | required? | meaning |
|--------|-----------|---------|
| `sample` | yes | unique sample name. |
| `fastq_1` | yes | read 1 / single-end fastq(.gz). |
| `fastq_2` | no | read 2 fastq(.gz); leave empty for single-end. Paired and single-end samples can be mixed. |
| `strandedness` | no | per-sample override: `auto`, `forward`, `reverse`, `unstranded`, or TEtranscripts' `no`. Blank = `auto`. |
| `condition` | no | biological group label. If present, TEtranscripts differential analysis runs automatically; if absent, it's skipped. |

## Test profile

`config/test.yaml` runs the whole workflow end-to-end against a bundled
synthetic dataset in `.tests/` (10kb genome, one gene, one TE, 4 paired-end
samples). Useful for confirming your setup or smoke-testing rule edits:

```bash
snakemake --configfile config/test.yaml --cores 4 -n   # dry-run
snakemake --configfile config/test.yaml --cores 4      # run
```

## Useful partial targets

```bash
snakemake --cores 8 star_index_only      # just build the STAR index
snakemake --cores 8 strandedness_only    # align + auto-detect strandedness only
```

## HPC / SLURM

Rule-level threads/memory/runtime defaults live in
`workflow/default-config/resources.yaml`. Override any of them by creating
`config/resources.yaml` with just the keys you want to change (partial
overrides are deep-merged). Unlisted rules fall back to conservative defaults
(1 thread / 4 GB / 60 min).

To submit to SLURM:

```bash
pip install snakemake-executor-plugin-slurm
snakemake --workflow-profile profiles/slurm
```

Edit `slurm_partition` / `slurm_account` in `profiles/slurm/config.yaml` for
your cluster. Per-invocation overrides are also possible, e.g.
`--set-resources star_align:mem_mb=64000`. Every tool runs from the shared
conda env, so make sure it's visible from the compute nodes (conda envs are
self-contained; if nodes can't read your home dir, install it on shared
storage with `conda env create -p /shared/path/env` and adjust your `PATH`).

## Tool versions

Every tool version is pinned independently with defaults shipped in
`workflow/default-config/versions.yaml` (STAR 2.7.11b, samtools 1.22.1, RSeQC
5.0.4, MultiQC 1.33, TEtranscripts 2.2.4, DESeq2 1.46.0, UCSC tools 482).
Override any by creating `config/versions.yaml` with just the keys you want to
change:

```yaml
versions:
  star: "2.7.10b"
```

The tools are installed into the shared conda env (`environment.yaml`) at env
creation, with the same pins as `versions.yaml` — keep the two files in sync
when bumping a version, then `conda env update -f environment.yaml` to apply.
`samtools` is pinned to 1.22.1 rather than newer: STAR 2.7.11b (the final STAR
release) links against `htslib <1.23`, which is incompatible with samtools
1.23+ in a single environment.

As a fallback, `--sdm conda` is still supported: Snakemake then generates one
small per-tool env per rule from `versions.yaml` at parse time instead of using
the shared env. This workflow does not use snakemake-wrappers, since a wrapper
tag pins all its tools' versions together.

One exception to the conda scheme: **TEtranscripts is installed from PyPI**
(`TEtranscripts==2.2.4`), not bioconda. The bioconda recipe's runtime deps pin
an ancient `bioconductor-deseq` (DESeq v1) that can only coexist with R 4.0-era
packages, making the conda package unsolvable alongside a modern DESeq2/R on
any platform. TEtranscripts only uses DESeq2 at runtime, so the environment
conda-installs `deseq2`/`r_base` (pins still apply) and pip-installs the
pure-Python TEtranscripts package.

## Without conda

If STAR/samtools/RSeQC/MultiQC/TEtranscripts are already on your `PATH` (e.g.
via `module load` on a cluster), just run `snakemake` without the shared env
activated — Snakemake runs the plain shell commands. You then own the version
matching; the shared conda env exists to guarantee it automatically.

## How strandedness is resolved

1. Per-sample value from the `strandedness` column, if filled in; otherwise `auto`.
2. Only `auto` samples go through RSeQC (`samtools sort/index` -> GTF converted
   to BED12 -> `infer_experiment.py` -> `determine_strandedness.py`). The
   forward- or reverse-strand fraction must reach `strandedness.min_fraction`
   to call `forward`/`reverse`; otherwise the library is called `no`
   (unstranded).
3. For `tetranscripts_diffexp`, all samples in a contrast must resolve to the
   same strandedness value or the run fails with a clear error.

TEcount always receives raw, unsorted STAR output, per the TEtranscripts
authors' recommendation.

## Automatic differential analysis

If `samples.csv` has a `condition` column, the workflow builds one
TEtranscripts contrast for every pairwise combination of distinct condition
values, named `{later}_vs_{earlier}` alphabetically (e.g. `control`/`treatment`
-> `treatment_vs_control`). If the column is absent, the workflow stops after
per-sample TEcount quantification.

## Output layout

```
results/
├── star/{sample}/Aligned.out.bam                       # unsorted, fed to TEcount
├── star/{sample}/Aligned.sortedByCoord.out.bam(.bai)   # for RSeQC/QC
├── rseqc/{sample}/infer_experiment.txt                 # raw RSeQC report
├── rseqc/{sample}/strandedness.txt                     # resolved no/forward/reverse
├── tecount/{sample}/{sample}.cntTable                  # per-sample gene+TE counts
├── tetranscripts/{contrast}/                           # only if "condition" column present
│   ├── {contrast}.cntTable
│   ├── {contrast}_DESeq2.R
│   ├── {contrast}_gene_TE_analysis.txt                 # full DESeq2 results
│   └── {contrast}_sigdiff_gene_TE.txt                  # significant hits only
└── qc/multiqc_report.html
```

Per-rule logs mirror this structure under `logs/` (plus
`logs/config_resolution.log`, which records how `sjdb_overhang` and gzipped
reference files were resolved).

## Notes

- Gzipped reference files (`.fa.gz` / `.gtf.gz`) are decompressed automatically;
  gzipped fastq files are read natively by STAR. Mix and match freely.
- TEtranscripts/DESeq2 needs at least 2 replicates per group in a contrast.
- STAR indexing and TEtranscripts are memory-hungry (TEtranscripts: ~20-30 GB
  recommended for human data).
