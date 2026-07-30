process INFER_STRANDEDNESS {
  label 'rseqc'
  publishDir "${params.outdir}/rseqc/${id}", mode: 'copy', pattern: "infer_experiment.txt"

  input:
  tuple val(id), path(sorted_bam), path(bai), val(condition), val(stranded_mode), path(bed12)

  output:
  tuple val(id), path(sorted_bam), path("infer_experiment.txt"), val(stranded_mode), val(condition), emit: inferred

  script:
  """
  infer_experiment.py \\
    --input-file $sorted_bam \\
    --refgene $bed12 \\
    > infer_experiment.txt 2> infer_experiment.log
  """

  stub:
  """
  echo -e 'Fraction of reads explained by "++,--": 0.85\\nFraction of reads explained by "+-,-+": 0.15' > infer_experiment.txt
  """
}
