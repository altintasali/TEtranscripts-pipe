# rnaseq-star-tetranscripts

[github.com/altintasali/rnaseq-star-tetranscripts](https://github.com/altintasali/rnaseq-star-tetranscripts)

A Snakemake workflow that quantifies genes **and** transposable elements from
RNA-seq data with [TEtranscripts/TEcount](https://github.com/mhammell-laboratory/TEtranscripts),
starting from STAR alignment and auto-detecting library strandedness with
RSeQC's `infer_experiment.py`.

```
fastq --> STAR align --> [samtools sort/index --> RSeQC infer_experiment --> strandedness] --> TEcount (per sample)
                                                                                            --> TEtranscripts (per contrast, DESeq2, only if "condition" column present)
                                                                             --> MultiQC
```

## Quick start

```
# 0. Get the workflow
git clone https://github.com/altintasali/rnaseq-star-tetranscripts.git
cd rnaseq-star-tetranscripts

# 1. Install Snakemake + this workflow's Python dependencies (see "Setup"
#    below for macOS/Homebrew/HPC-specific instructions if this doesn't fit)
mamba env create -f environment.yaml
conda activate rnaseq-star-tetranscripts

# 2. Build the container image with all bioinformatics tools (STAR, samtools,
#    TEtranscripts, RSeQC, MultiQC, UCSC tools). Requires Docker installed.
docker build -t rnaseq-star-tetranscripts .

# 3. Try it against the bundled tiny synthetic dataset -- no real genome/
#    reads needed, runs end-to-end in a couple of minutes (see "Test profile")
snakemake --configfile config/test.yaml --sdm docker --cores 4

# 4. Run on your own data: edit config/config.yaml (genome fasta/gtf/te_gtf
#    paths) and config/samples.csv (your samples), then:
snakemake --sdm docker --cores 16
```

That's the whole workflow. Everything else in this README is reference
material for when you need to go beyond the defaults: tool versions,
per-sample strandedness overrides, differential analysis, HPC/SLURM, etc.

Configuration follows nf-core's "sensible defaults, minimal user-facing
config" philosophy. `workflow/default-config/{versions,resources}.yaml`
ship built-in defaults as part of the workflow itself -- you never need to
open them for a normal run. The only files you touch to run an analysis
are `config/config.yaml` (ref paths, STAR/TEtranscripts settings) and
`config/samples.csv` (the design file). If you ever do need to change a
tool version or a rule's CPU/memory, create `config/versions.yaml` and/or
`config/resources.yaml` yourself with just the keys you want to change --
Snakemake deep-merges them on top of the built-in defaults, so partial
overrides are fine. See "Tool versions" and "HPC / SLURM" below.

**This branch replaces per-rule conda environments with a single container
image (Docker / Singularity).** Every rule runs inside the container
defined by `Dockerfile` at the repo root. On local machines / CI, use
`--sdm docker`; on HPC clusters with Singularity/Apptainer, push the
image to a registry and set `--sdm apptainer` (see "HPC / SLURM" below).

## Setup

1. Install Snakemake **and its Python dependencies** into one environment.
   Pick the section below that matches your situation.

   ### All platforms: conda/mamba (for Snakemake + Python deps only)
   The only purpose of this environment is to provide Snakemake and the
   Python packages it needs at workflow-parse time (`pandas`, `pyyaml`,
   `jsonschema`). All bioinformatics tools (STAR, samtools, TEtranscripts,
   etc.) live in a **container image** (see `Dockerfile`) and are run via
   `--sdm docker` or `--sdm apptainer` -- conda is no longer used for
   per-rule tool environments.

   ```
   mamba env create -f environment.yaml
   conda activate rnaseq-star-tetranscripts
   docker build -t rnaseq-star-tetranscripts .
   ```
   No conda/mamba yet? Install [Miniforge](https://github.com/conda-forge/miniforge)
   (conda + mamba in one installer) -- on macOS that's also `brew install
   miniforge` followed by `conda init "$(basename "${SHELL}")"` and
   restarting your terminal.

   ### macOS via Homebrew (`brew install snakemake`)
   Homebrew's `snakemake` formula installs into its own isolated
   environment. This works fine -- just add the missing dependencies:
   ```
   $(brew --prefix snakemake)/libexec/bin/pip install pandas pyyaml jsonschema
   ```
   Then build the container image:
   ```
   docker build -t rnaseq-star-tetranscripts .
   ```

    ### HPC clusters / Singularity + SLURM
    Most HPC clusters provide Singularity/Apptainer but not Docker. The
    workflow runs identically under Apptainer:
    1. On a machine with Docker (your laptop, a login node, CI), build:
       ```
       docker build -t rnaseq-star-tetranscripts .
       ```
    2. Push the image to a container registry:
       ```
       docker tag rnaseq-star-tetranscripts your-user/rnaseq-star-tetranscripts:latest
       docker push your-user/rnaseq-star-tetranscripts:latest
       ```
    3. Update `CONTAINER_IMAGE` in `workflow/rules/common.smk` to point at
       your pushed image:
       ```python
       CONTAINER_IMAGE = "docker://your-user/rnaseq-star-tetranscripts:latest"
       ```
       A pre-built image is available on Docker Hub:
       ```
       CONTAINER_IMAGE = "docker://altintasali/rnaseq-star-tetranscripts:latest"
       ```
    4. Load Snakemake (via `module load` or your own conda env) and run:
       ```
       module load miniconda3 apptainer
       snakemake --workflow-profile profiles/slurm --sdm apptainer --cores 32
       ```
   A few HPC-specific things worth knowing:
   - Compute nodes on most clusters cannot reach the internet to pull
     container images -- pull the image on a **login node** first so the
     cached SIF file is available to compute nodes.
   - If your cluster shares a filesystem across nodes, Apptainer's cache
     (`~/.apptainer/cache`) is readable from any compute node; otherwise
     set `APPTAINER_CACHEDIR` to a shared path.

2. Get the workflow:
   ```
   git clone https://github.com/altintasali/rnaseq-star-tetranscripts.git
   cd rnaseq-star-tetranscripts
   ```
   Or deploy it into an existing project via
   [`snakedeploy`](https://snakedeploy.readthedocs.io/):
   ```
   snakedeploy deploy-workflow https://github.com/altintasali/rnaseq-star-tetranscripts . --tag main
   ```
3. Provide reference files and edit `config/config.yaml` (the minimum you
   need to touch to run an analysis):
   - `ref.fasta` / `ref.gtf`: genome FASTA and gene GTF (Ensembl/GENCODE).
     Plain or gzipped (`.fa.gz`/`.gtf.gz`) both work -- gzipped files are
     decompressed automatically.
   - `ref.te_gtf`: the **curated** TE GTF from the TEtranscripts authors --
     download the file matching your genome build from
     https://www.mghlab.org/software/tetranscripts (a generic RepeatMasker
     GTF will *not* work correctly with TEtranscripts). Plain or gzipped
     also both work here.
   - `ref.sjdb_overhang`: leave as `auto` (default) to detect it from the
     sample sheet's fastq read lengths automatically, or set an explicit
     integer to pin it. See "STAR sjdbOverhang auto-detection" below.
   - `star.index`: where to (re)build the STAR index.
   - `strandedness.min_fraction`: confidence threshold for RSeQC
     auto-detection (per-sample strandedness itself is set in
     `samples.csv`, not here -- see below).
   - `tetranscripts.*`: TEcount/TEtranscripts options (mode, padj, foldchange...).

   `config/versions.yaml` and `config/resources.yaml` don't exist by
   default and you don't need to create them -- built-in defaults for both
   already ship in `workflow/default-config/`. Create either file yourself,
   with just the keys you want to change, only if you need to pin a
   different tool version ("Tool versions" below) or retune CPU/memory for
   your cluster ("HPC / SLURM" below).
4. Fill in `config/samples.csv` (the design file), columns **in this exact
   order** -- this matches the nf-core/rnaseq samplesheet layout, so an
   nf-core samplesheet can be used directly:

   | column       | required? | meaning                                                        |
   |--------------|-----------|------------------------------------------------------------------|
   | sample       | yes       | unique sample name, used throughout the workflow                  |
   | fastq_1      | yes       | path to read 1 / single-end fastq(.gz)                            |
   | fastq_2      | no        | path to read 2 fastq(.gz). Leave empty for a single-end sample -- paired- and single-end samples can be freely mixed in one sheet. |
   | strandedness | no        | per-sample override: `auto`, `forward`, `reverse`, or `unstranded` (nf-core vocabulary; TEtranscripts' own `no` is also accepted). Leave empty to default to `auto` for that sample. Mixing values across samples is fine. |
   | condition    | no        | biological group label. **If this column is present**, TEtranscripts differential analysis runs automatically (see below). **If absent, differential analysis is skipped** and the workflow stops after per-sample TEcount quantification. |

## Run

```
# sanity check
snakemake -n --sdm docker

# full run
snakemake --sdm docker --cores 16
```

Useful partial targets:
```
snakemake --sdm docker --cores 8 star_index_only      # just build the STAR index
snakemake --sdm docker --cores 8 strandedness_only    # align + auto-detect strandedness only
```

On HPC with Singularity/Apptainer, replace `--sdm docker` with `--sdm apptainer`
(or `--sdm singularity` for older versions).

## Test profile

`config/test.yaml` runs the whole workflow (STAR index/align, RSeQC
auto-detection, TEcount, TEtranscripts diffexp, MultiQC) end-to-end against
a tiny bundled synthetic dataset in `.tests/` -- a 10kb single-contig
genome with one gene and one TE, 4 paired-end samples (2
treatment/2 control), all generated deterministically by
`.tests/generate_test_data.py`. No real genome/reads needed; useful for
confirming your Docker/container setup works, or as a smoke test after
editing any rule:

```
snakemake --configfile config/test.yaml --sdm docker --cores 4 -n   # dry-run
snakemake --configfile config/test.yaml --sdm docker --cores 4      # run
```

`config/test.yaml` is just another config layer -- passed on the command
line, it's deep-merged on top of everything else (including your own
`config/config.yaml`), so it only needs to override the handful of keys
that differ for the toy dataset: `samples`/`ref.*` point at `.tests/`,
`star.index_extra` sets `--genomeSAindexNbases 5` (STAR requires shrinking
this for genomes well under its default ~2^30bp assumption), and
`resources` shrinks every rule's CPU/memory for a quick laptop/CI run. To
regenerate or resize the test data, edit and rerun
`.tests/generate_test_data.py`.

## HPC / SLURM

Every rule's threads, memory, and runtime come from
`workflow/default-config/resources.yaml` (one entry per rule; see that
file for the current defaults and per-rule notes on why they're sized that
way), via a small `get_resources()` lookup in `common.smk`. To retune
CPU/memory for your cluster or dataset size, create `config/resources.yaml`
with just the rules/keys you want to change (it overrides the built-in
defaults; see "Configuration" above) -- no need to touch the rule files
themselves. A rule left out of both files (or missing a key) falls back to
a conservative `{threads: 1, mem_mb: 4000, runtime: 60}` default rather
than failing.

To submit to a SLURM cluster, install the executor plugin, build/push the
container image, and update `CONTAINER_IMAGE` in
`workflow/rules/common.smk` to point at your registry (see "Setup > HPC
clusters" above), then use the included workflow profile:

```
pip install snakemake-executor-plugin-slurm

snakemake --workflow-profile profiles/slurm --sdm apptainer
```

`profiles/slurm/config.yaml` sets `executor: slurm`, a job cap (`jobs: 50`),
and fallback `default-resources` (partition/account/mem_mb/runtime) for
anything not covered by `resources.yaml`. Edit `slurm_partition` and
`slurm_account` in that file for your cluster before running. Each
Snakemake job becomes one `sbatch` submission sized from its rule's
`threads`/`mem_mb`/`runtime` -- e.g. `star_align` submits as a 12-CPU,
48GB, 4-hour job per sample, `tetranscripts_diffexp` as a 4-CPU, 32GB,
6-hour job per contrast, and so on.

You can override any of this per-invocation without editing files, e.g. to
bump memory for one run:
```
snakemake --workflow-profile profiles/slurm --sdm apptainer \
  --set-resources star_align:mem_mb=64000
```

Without SLURM, everything above still applies locally -- `--cores N` just
caps how many rules run in parallel on the machine you're on, respecting
each rule's `threads`.

## Tool versions -- specify every tool independently, no manual installs

Every tool version is pinned individually, with defaults shipped in
`workflow/default-config/versions.yaml` (part of the workflow, not meant
to be edited there):

```yaml
versions:
  star: "2.7.11b"
  samtools: "1.23.1"
  rseqc: "5.0.4"
  multiqc: "1.33"
  tetranscripts: "2.2.4"
  ucsc_gtftogenepred: "482"
  ucsc_genepredtobed: "482"
```

To change one or more, create `config/versions.yaml` with just the keys
you want to override, e.g.:

```yaml
versions:
  star: "2.7.10b"
```

`workflow/Snakefile` loads `workflow/default-config/versions.yaml` first,
then `config/config.yaml`, then `config/versions.yaml` only if it exists;
Snakemake deep-merges all of them into one `config` dict, so a partial
override is fine and the rest of the workflow doesn't need to know they're
separate files.

The `Dockerfile` at the repo root reads these version strings at build
time and installs each tool at its pinned version. To change a version,
edit `config/versions.yaml` (or the defaults in
`workflow/default-config/versions.yaml`), then rebuild the container
image:

```
docker build --no-cache -t rnaseq-star-tetranscripts .
```

**To match a specific external pipeline's tool versions** (e.g. a particular
nf-core/rnaseq release), pull the exact numbers from that run's
`pipeline_info/software_versions.yml`, or the "Software Versions" section of
its MultiQC report, and paste them into `config/versions.yaml`. The
shipped defaults are already set to match a specific nf-core/rnaseq run;
override them if you're targeting a different one.

## Container-based execution

This workflow uses a **single** container image (defined by `Dockerfile`
at the repo root) that bundles every tool: STAR, samtools, RSeQC,
MultiQC, TEtranscripts/TEcount, UCSC tools, and R/Bioconductor. Every
rule's `shell:` command runs inside this image via:

- `--sdm docker` on local machines / CI (Docker available)
- `--sdm apptainer` on HPC clusters (Apptainer/Singularity available)

The `CONTAINER_IMAGE` variable in `workflow/rules/common.smk` controls
which image is used. By default it points to `dockerfile://Dockerfile`
(build from the local Dockerfile). For production/deployment, build the
image, push it to a container registry, and update `CONTAINER_IMAGE` to
the registry URI (e.g. `docker://your-user/image:tag`). A pre-built image
is available on Docker Hub as `altintasali/rnaseq-star-tetranscripts:latest`.

Versions of every tool in the image are pinned in
`workflow/default-config/versions.yaml`. After changing a version, rebuild
the image:

```
docker build --no-cache -t rnaseq-star-tetranscripts .
```

### Running without a container manager

If neither Docker nor Apptainer/Singularity is available (e.g. testing on
a bare-metal server), drop `--sdm docker`/`--sdm apptainer` and Snakemake
will run the plain shell commands against whatever is on your `PATH` --
you are responsible for matching the tool versions to those in
`versions.yaml`. The `container:` directive in the workflow is silently
ignored when no `--sdm` flag is passed.

## How strandedness is resolved

For each sample, the effective strandedness mode is resolved in this order:
1. The sample sheet's `strandedness` column for that sample, if filled in
   (`auto`/`forward`/`reverse`/`unstranded`, or TEtranscripts' own `no`).
2. Otherwise, `auto` by default.

Only samples whose *effective* mode ends up `auto` go through RSeQC
auto-detection; the rest use their fixed value directly, so a sample sheet
can freely mix explicit and auto-detected samples. For samples that need it:
1. STAR's BAM is sorted + indexed (`samtools sort`/`index`).
2. `ref.gtf` is converted to a BED12 gene model once (`gtfToGenePred` +
   `genePredToBed`, UCSC tools) for RSeQC.
3. `infer_experiment.py` runs on the sorted BAM against that BED12 model.
4. `workflow/scripts/determine_strandedness.py` parses the two "Fraction of
   reads explained by ..." lines and maps them to TEtranscripts/TEcount's
   `--stranded` values:
   - forward-type fraction dominates (>= `strandedness.min_fraction`) -> `forward`
   - reverse-type fraction dominates (>= `strandedness.min_fraction`) -> `reverse`
   - neither dominates -> `no` (unstranded)

For the `tetranscripts_diffexp` rule (which pools several BAMs into one
DESeq2 run), the workflow requires every sample in a contrast to resolve to
the *same* final strandedness value (whether fixed or auto-detected), and
fails with a clear error otherwise -- that usually means samples in the
contrast were prepped with different library kits and shouldn't be pooled
as-is.

TEcount itself always receives raw, unsorted STAR output
(`Aligned.out.bam`), per the TEtranscripts authors' recommendation
(unsorted or queryname-sorted input only).

## Automatic differential analysis

If `config/samples.csv` has a `condition` column filled in for every sample,
the workflow builds one TEtranscripts contrast for **every pairwise
combination** of distinct condition values, named `{later}_vs_{earlier}` in
alphabetical order (e.g. `control`/`treatment` -> one contrast,
`treatment_vs_control`; `control`/`kd1`/`kd2` -> three contrasts:
`kd1_vs_control`, `kd2_vs_control`, `kd2_vs_kd1`). There is nothing to
configure -- just add the column. If you want a specific pair labeled
"treatment" vs. "control" rather than alphabetical order, name your
condition values accordingly (e.g. `a_control` / `b_treatment`), or edit
`_build_contrasts()` in `workflow/rules/common.smk` for full manual control.

If the `condition` column is absent, `CONTRASTS` is empty and the workflow
stops after per-sample `TEcount` quantification -- no TEtranscripts/DESeq2
step is requested or run.

## STAR sjdbOverhang auto-detection

`ref.sjdb_overhang` defaults to `auto`. At workflow parse time (before any
rule runs), `common.smk` peeks at the first ~20 reads of every `fastq_1`/
`fastq_2` file in the sample sheet (gzip-aware) and sets
`sjdbOverhang = max(read length) - 1` across all of them, per STAR's own
recommendation -- one value is used for the whole shared genome index. This
means the fastq files need to actually exist/be reachable at parse time; if
none can be read (e.g. paths are placeholders, or reads aren't downloaded
yet), the workflow fails immediately with a clear error rather than
building an index with a guessed value. Set `ref.sjdb_overhang` to an
explicit integer in `config.yaml` to skip auto-detection entirely (e.g. if
you want to build the index before fastq files are in place).

The resolved value (and its source) is written to
`logs/config_resolution.log` on every run for a quick sanity check.

## Gzipped reference files

`ref.fasta`, `ref.gtf`, and `ref.te_gtf` can each be plain text or gzipped
(`.fa.gz`/`.gtf.gz`) independently -- mix and match freely. Any of them
ending in `.gz` is decompressed once by the `gunzip_reference` rule into
`resources/decompressed/`, since STAR's `--genomeFastaFiles`/
`--sjdbGTFfile` and TEcount/TEtranscripts' `--GTF`/`--TE` all expect plain
text. Every other rule always refers to the resolved (guaranteed
uncompressed) path, so nothing downstream needs to know or care whether the
original file was gzipped. `logs/config_resolution.log` lists which files
(if any) triggered decompression on the current run.

fastq files are unaffected by this -- STAR reads gzipped fastq natively
(the workflow auto-adds `--readFilesCommand zcat` when `fastq_1` ends in
`.gz`), so paired/single-end and compressed/uncompressed reads can all be
mixed freely in one sample sheet.

## Logs

Every step that runs an external tool writes its own log file under
`logs/`, mirroring the rule/sample structure of `results/`, e.g.:

```
logs/
├── config_resolution.log             # sjdb_overhang + gzip decompression decisions (see above)
├── gunzip/{stem}.log                 # only for reference files that were gzipped
├── star/index.log
├── star/align/{sample}.log
├── samtools/sort/{sample}.log
├── samtools/index/{sample}.log
├── rseqc/gtf_to_genepred.log
├── rseqc/genepred_to_bed12.log
├── rseqc/infer_experiment/{sample}.log
├── rseqc/determine_strandedness/{sample}.log
├── tecount/{sample}.log
├── tetranscripts/{contrast}.log
└── multiqc.log
```

These are plain stdout/stderr redirects from each tool (or, for the Python
scripts, `snakemake.log`), so they're the first place to look when a rule
fails or a step's `results/` output looks wrong.

## Output layout

```
results/
├── star/{sample}/Aligned.out.bam                       # unsorted, fed to TEcount/TEtranscripts
├── star/{sample}/Aligned.sortedByCoord.out.bam(.bai)    # for RSeQC/QC only
├── rseqc/{sample}/infer_experiment.txt                  # raw RSeQC report
├── rseqc/{sample}/strandedness.txt                      # resolved no/forward/reverse
├── tecount/{sample}/{sample}.cntTable                   # per-sample gene+TE counts
├── tetranscripts/{contrast}/                            # only if "condition" column present
│   ├── {contrast}.cntTable
│   ├── {contrast}_DESeq2.R
│   ├── {contrast}_gene_TE_analysis.txt                  # full DESeq2 results
│   └── {contrast}_sigdiff_gene_TE.txt                   # significant hits only
└── qc/multiqc_report.html
```

## Notes / things worth reviewing before a production run

- `star.extra` in the config bakes in the multi-mapper settings the
  TEtranscripts authors recommend (`--outFilterMultimapNmax 100
  --winAnchorMultimapNmax 100`) -- tune these for your organism/experiment.
- The STAR index build and TEcount/TEtranscripts steps are memory-hungry
  (STAR indexing: genome-size dependent; TEtranscripts: ~20-30GB recommended
  for human at 20-30M mapped reads). Add `resources: mem_mb=...` to the
  relevant rules and use a cluster profile for HPC/cloud execution.
- `TEtranscripts`/DESeq2 needs at least 2 replicates per group in a contrast.
