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
