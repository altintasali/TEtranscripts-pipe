# TO-DO List

- `tetranscripts-pipe` script (`init` and `samples`) should be more flexible for creating the samplesheet.csv. It should be written to current directory
- make all "write_igv_to_bed" rules "true"
- strandedness check MultiQC text is wrong when strandedness is declared in samples.csv
  (`workflow/scripts/strandedness_check_mqc.py:135-143`). For every "undetermined" sample it
  says `Quantification ran with "no"`, but that's only true when the declared mode is `auto`.
  When the sheet sets `forward`/`reverse`, the declared value is used (`get_strandedness_param`),
  and the table's own "used" column correctly shows it, so the text contradicts the table.
  Fix: split `n_undetermined` by declared mode (auto -> "ran with no"; declared -> "ran with the
  declared value; RSeQC alone could not confirm it") and add a test. Seen in
  20261111_Setdb1-KD-hESC (5/6 samples at 75-81% antisense, declared `reverse`).
- changing a sample's strandedness in samples.csv does not rerun finished strand-dependent jobs.
  Rules that take `--stranded` from `get_strandedness_param` (`tecount`, `telocal`,
  `chimera_chimeric_reads_classify`) keep outputs made with the old value, so results silently
  mix strand settings. Seen in 20261111_Setdb1-KD-hESC: sheet changed `auto` -> `reverse` and the
  run restarted; only jobs that hadn't finished were redone, and N_DMSO_R3/N_DOX_R2
  `chimera_chimeric_reads_classify` outputs stayed from the `no` run (needed `--forcerun`).
  Likely cause: the params trigger doesn't see the change because the value comes from a params
  function. Fix idea: write the resolved strandedness per sample to a file (e.g.
  `results/rseqc/{sample}_strandedness_used.txt`, rewritten only when it changes) and make it an
  input of these rules, so a change reruns them via mtime; add a test.
- MultiQC report collapses replicates whose sample names contain `_R1`/`_R2`
  (`workflow/default-config/multiqc_config.yaml:56-61`). `extra_fn_clean_exts` lists bare
  `"_R1"`/`"_R2"` (meant to merge per-read FastQC/cutadapt reports), but bare entries are MultiQC
  `truncate` rules: the name is cut at the *first* match anywhere. So `N_DMSO_R1` and `N_DMSO_R2`
  both become `N_DMSO` (one overwrites the other), while `N_DMSO_R3_R1` -> `N_DMSO_R3` survives;
  the report showed 4 of 6 samples in every section, including the pipeline's own PCA custom
  content. Counts/outputs are unaffected. Seen in 20261111_Setdb1-KD-hESC (worked around by
  renaming samples to `rep1`-`rep3`); nf-core/rnaseq is fine with the same names because it
  suffixes reads `_1`/`_2` and never cleans `_R1`/`_R2`.
  An end-anchored `type: regex, pattern: "_R[12]$"` is NOT enough: sample-level outputs
  (STAR logs, custom content) are named by the bare sample, so `N_DMSO_R1` would still lose `_R1`.
  Fix: drop `_R1`/`_R2` from `extra_fn_clean_exts` and let the `multiqc` rule (`rules/qc.smk`)
  write a `--replace-names` TSV from the sample sheet, mapping each per-read name to its sample
  (`{sample}_R1` -> `{sample}`, `{sample}_R2` -> `{sample}`; or "`{sample} Read 1/2`" with
  `table_sample_merge`, as nf-core does), with `sample_names_replace_exact: true`. Exact mapping
  from known sample names cannot collide with a replicate tag. Add a test with samples named
  `X_R1`, `X_R2`, `X_R3` asserting all three appear in `multiqc_general_stats.txt`.
