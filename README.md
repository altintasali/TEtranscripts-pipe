# rna-seq-star-tetranscripts (custom)

A Snakemake workflow that quantifies genes **and** transposable elements from
RNA-seq data with [TEtranscripts/TEcount](https://github.com/mhammell-laboratory/TEtranscripts),
starting from STAR alignment and auto-detecting library strandedness with
RSeQC's `infer_experiment.py`.

```
fastq --> STAR align --> [samtools sort/index --> RSeQC infer_experiment --> strandedness] --> TEcount (per sample)
                                                                                            --> TEtranscripts (per contrast, DESeq2, only if "condition" column present)
                                                                             --> MultiQC
```

Configuration is split across three files, all loaded and deep-merged by
`workflow/Snakefile`: `config/config.yaml` for day-to-day analysis settings
(the one you actually edit per run), `config/versions.yaml` for tool/package
versions (rarely needs touching), and `config/resources.yaml` for per-rule
CPU/memory/runtime (only relevant for cluster/HPC tuning -- see "HPC /
SLURM" below). Every tool version (STAR, samtools, RSeQC, MultiQC,
TEtranscripts/TEcount, DESeq2, UCSC gtfToGenePred/genePredToBed) is pinned
independently and rendered into its own conda environment automatically --
see "Tool versions" below.

## Setup

1. Install Snakemake (>=8) via conda/mamba:
   ```
   mamba create -c conda-forge -c bioconda -n snakemake snakemake
   conda activate snakemake
   ```
2. Get the workflow (clone this directory, or `snakedeploy` it if you push it
   to your own GitHub repo).
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

   `config/versions.yaml` holds tool/package versions and rarely needs
   editing -- see "Tool versions" below. `config/resources.yaml` holds
   per-rule CPU/memory/runtime and only matters for cluster/HPC tuning --
   see "HPC / SLURM" below.
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
snakemake -n --sdm conda

# full run
snakemake --sdm conda --cores 16
```

Useful partial targets:
```
snakemake --sdm conda --cores 8 star_index_only      # just build the STAR index
snakemake --sdm conda --cores 8 strandedness_only    # align + auto-detect strandedness only
```

## HPC / SLURM

Every rule's threads, memory, and runtime come from `config/resources.yaml`
(one entry per rule; see that file for the current defaults and per-rule
notes on why they're sized that way), via a small `get_resources()` lookup
in `common.smk` -- edit that one file to retune CPU/memory for your cluster
or dataset size, no need to touch the rule files themselves. A rule left
out of `resources.yaml` (or missing a key) falls back to a conservative
`{threads: 1, mem_mb: 4000, runtime: 60}` default rather than failing.

To submit to a SLURM cluster, install the executor plugin and use the
included workflow profile:

```
pip install snakemake-executor-plugin-slurm

snakemake --workflow-profile profiles/slurm --sdm conda
```

`profiles/slurm/config.yaml` sets `executor: slurm`, a job cap (`jobs: 50`),
and fallback `default-resources` (partition/account/mem_mb/runtime) for
anything not covered by `resources.yaml`. Edit `slurm_partition` and
`slurm_account` in that file for your cluster before running. Each
Snakemake job becomes one `sbatch` submission sized from its rule's
`threads`/`mem_mb`/`runtime` (`config/resources.yaml`) -- e.g. `star_align`
submits as a 12-CPU, 48GB, 4-hour job per sample, `tetranscripts_diffexp`
as a 4-CPU, 32GB, 6-hour job per contrast, and so on.

You can override any of this per-invocation without editing files, e.g. to
bump memory for one run:
```
snakemake --workflow-profile profiles/slurm --sdm conda \
  --set-resources star_align:mem_mb=64000
```

Without SLURM, everything above still applies locally -- `--cores N` just
caps how many rules run in parallel on the machine you're on, respecting
each rule's `threads`.

## Tool versions -- specify every tool independently, no manual installs

Every tool version is pinned individually in `config/versions.yaml`,
separate from the analysis settings in `config/config.yaml`:

```yaml
versions:
  star: "2.7.11b"
  samtools: "1.23.1"
  rseqc: "5.0.4"
  multiqc: "1.33"
  tetranscripts: "2.2.4"
  deseq2: "1.46.0"
  ucsc_gtftogenepred: "487"
  ucsc_genepredtobed: "487"
```

`workflow/Snakefile` loads all three files (`configfile: "config/config.yaml"`,
then `"config/versions.yaml"`, then `"config/resources.yaml"`); Snakemake
deep-merges them into
one `config` dict, so the rest of the workflow doesn't need to know they're
separate files.

This workflow does **not** use snakemake-wrappers for STAR/samtools/RSeQC/
MultiQC (or anything else) -- deliberately. A wrapper's git tag bundles
*all* of its tools' versions together, so you can't independently ask for
"STAR 2.7.11b + samtools 1.23.1" if the wrapper tag that has one doesn't
also have the other. Instead, `common.smk` renders each version string
above into its own tiny conda env file at parse time
(`workflow/envs/generated/{star,samtools,rseqc,multiqc,tetranscripts,ucsc_tools}.yaml`),
and every rule runs the plain command directly (`STAR ...`, `samtools sort
...`, `infer_experiment.py ...`, `multiqc ...`) in that env.

Running with `--sdm conda` (or the older `--use-conda`) makes Snakemake
create/download the exact matching environment automatically the first time
you run -- there's nothing to install by hand. Just edit a version string
and rerun.

**To match a specific external pipeline's tool versions** (e.g. a particular
nf-core/rnaseq release), pull the exact numbers from that run's
`pipeline_info/software_versions.yml`, or the "Software Versions" section of
its MultiQC report, and paste them into `config/versions.yaml` above -- that
file is the authoritative source, since modern nf-core modules can resolve
tool versions dynamically from their containers rather than a single
grep-able pin in the module source. The defaults above are already set to
match a specific nf-core/rnaseq run; adjust them if you're targeting a
different one.


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
