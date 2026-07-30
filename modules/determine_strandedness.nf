process DETERMINE_STRANDEDNESS {
  label 'strandedness'
  publishDir "${params.outdir}/rseqc/${id}", mode: 'copy', pattern: "strandedness.txt"

  input:
  tuple val(id), path(sorted_bam), path(infer_txt), val(stranded_mode), val(condition)
  path strandedness_script

  output:
  tuple val(id), path(sorted_bam), path("strandedness.txt"), val(condition), emit: stranded

  script:
  """
  python $strandedness_script $infer_txt strandedness.txt ${params.strandedness_min_fraction}
  """

  stub:
  """
  echo "forward" > strandedness.txt
  """
}
