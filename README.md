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

Requires [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba; on
macOS also `brew install miniforge`).

```bash
git clone https://github.com/altintasali/rnaseq-star-tetranscripts.git
cd rnaseq-star-tetranscripts
mamba env create -f environment.yaml        # snakemake + its python deps
conda activate rnaseq-star-tetranscripts

# Smoke test on the bundled tiny synthetic dataset (no real genome/reads needed):
snakemake --configfile config/test.yaml --sdm conda --cores 4

# Run on your own data (edit config/config.yaml + config/samples.csv first):
snakemake --sdm conda --cores 16
```

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
snakemake --configfile config/test.yaml --sdm conda --cores 4 -n   # dry-run
snakemake --configfile config/test.yaml --sdm conda --cores 4      # run
```

## Useful partial targets

```bash
snakemake --sdm conda --cores 8 star_index_only      # just build the STAR index
snakemake --sdm conda --cores 8 strandedness_only    # align + auto-detect strandedness only
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
snakemake --workflow-profile profiles/slurm --sdm conda
```

Edit `slurm_partition` / `slurm_account` in `profiles/slurm/config.yaml` for
your cluster. Per-invocation overrides are also possible, e.g.
`--set-resources star_align:mem_mb=64000`. Build the per-rule conda envs once
on a login node with internet access first, and point `--conda-prefix` at
shared storage if compute nodes can't read the default cache.

## Tool versions

Every tool version is pinned independently with defaults shipped in
`workflow/default-config/versions.yaml` (STAR 2.7.11b, samtools 1.23.1, RSeQC
5.0.4, MultiQC 1.33, TEtranscripts 2.2.4, DESeq2 1.46.0, UCSC tools 487).
Override any by creating `config/versions.yaml` with just the keys you want to
change:

```yaml
versions:
  star: "2.7.10b"
```

`--sdm conda` creates the matching per-tool conda envs automatically on first
run — nothing to install by hand. (Without `--sdm conda`, the version pins are
ignored and whatever tools are on your `PATH` get used — see "Without conda"
below. This workflow does not use snakemake-wrappers, since a wrapper tag
pins all its tools' versions together.)

## Without conda

If STAR/samtools/RSeQC/MultiQC/TEtranscripts are already on your `PATH` (e.g.
via `module load` on a cluster), drop `--sdm conda` and Snakemake will run the
plain shell commands. You then own the version matching — `--sdm conda` exists
to guarantee it automatically.

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
