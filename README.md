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

Every tool version (STAR, samtools, RSeQC, MultiQC, TEtranscripts/TEcount,
DESeq2, UCSC gtfToGenePred/genePredToBed) is pinned independently in
`config["versions"]` and rendered into its own conda environment
automatically -- see "Tool versions" below.

## Setup

1. Install Snakemake (>=8) via conda/mamba:
   ```
   mamba create -c conda-forge -c bioconda -n snakemake snakemake
   conda activate snakemake
   ```
2. Get the workflow (clone this directory, or `snakedeploy` it if you push it
   to your own GitHub repo).
3. Provide reference files and edit `config/config.yaml`:
   - `ref.fasta` / `ref.gtf`: genome FASTA and gene GTF (Ensembl/GENCODE).
   - `ref.te_gtf`: the **curated** TE GTF from the TEtranscripts authors --
     download the file matching your genome build from
     https://www.mghlab.org/software/tetranscripts (a generic RepeatMasker
     GTF will *not* work correctly with TEtranscripts).
   - `ref.sjdb_overhang`: read length - 1.
   - `star.index`: where to (re)build the STAR index.
   - `strandedness.mode`: leave as `auto` to auto-detect per sample with
     RSeQC, or set `no`/`forward`/`reverse` to force one value.
   - `tetranscripts.*`: TEcount/TEtranscripts options (mode, padj, foldchange...).
   - `versions.*`: tool/package versions (see "Tool versions" below).
4. Fill in `config/samples.csv` (the design file), columns **in this order**:

   | column    | required? | meaning                                                        |
   |-----------|-----------|------------------------------------------------------------------|
   | sample    | yes       | unique sample name, used throughout the workflow                  |
   | fastq_1   | yes       | path to read 1 / single-end fastq(.gz)                            |
   | fastq_2   | no        | path to read 2 fastq(.gz). Leave empty for a single-end sample -- paired- and single-end samples can be freely mixed in one sheet. |
   | condition | no        | biological group label. **If this column is present**, TEtranscripts differential analysis runs automatically (see below). **If absent, differential analysis is skipped** and the workflow stops after per-sample TEcount quantification. |

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

## Tool versions -- specify every tool independently, no manual installs

Every tool version is pinned individually in `config["versions"]`:

```yaml
versions:
  star: "2.7.11b"
  samtools: "1.21"
  rseqc: "5.0.4"
  multiqc: "1.25.1"
  tetranscripts: "2.2.4"
  deseq2: "1.44.0"
  ucsc_gtftogenepred: "447"
  ucsc_genepredtobed: "447"
```

This workflow does **not** use snakemake-wrappers for STAR/samtools/RSeQC/
MultiQC (or anything else) -- deliberately. A wrapper's git tag bundles
*all* of its tools' versions together, so you can't independently ask for
"STAR 2.7.11b + samtools 1.21" if the wrapper tag that has one doesn't also
have the other. Instead, `common.smk` renders each version string above into
its own tiny conda env file at parse time
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
its MultiQC report, and paste them into `config["versions"]` above -- that
file is the authoritative source, since modern nf-core modules can resolve
tool versions dynamically from their containers rather than a single
grep-able pin in the module source. The defaults above are reasonable
current versions as a starting point, not a guaranteed match to any
particular pipeline release.

## How strandedness auto-detection works

For each sample, after STAR alignment the workflow:
1. Sorts + indexes the BAM (`samtools sort`/`index`).
2. Converts `ref.gtf` to a BED12 gene model once (`gtfToGenePred` +
   `genePredToBed`, UCSC tools) for RSeQC.
3. Runs `infer_experiment.py` (`bio/rseqc/infer_experiment` wrapper) on the
   sorted BAM against that BED12 model.
4. `workflow/scripts/determine_strandedness.py` parses the two "Fraction of
   reads explained by ..." lines and maps them to TEtranscripts/TEcount's
   `--stranded` values:
   - forward-type fraction dominates (>= `strandedness.min_fraction`) -> `forward`
   - reverse-type fraction dominates (>= `strandedness.min_fraction`) -> `reverse`
   - neither dominates -> `no` (unstranded)

For the `tetranscripts_diffexp` rule (which pools several BAMs into one
DESeq2 run), the workflow requires every sample in a contrast to have been
called with the *same* strandedness, and fails with a clear error otherwise
-- that usually means samples in the contrast were prepped with different
library kits and shouldn't be pooled as-is.

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
