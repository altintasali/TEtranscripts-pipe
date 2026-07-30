process GENEPRED_TO_BED12 {
  label 'genepred'

  input:
  path genepred

  output:
  path("annotation.bed12"), emit: bed12

  script:
  """
  genePredToBed $genepred annotation.bed12 > genepred_to_bed12.log 2>&1
  """

  stub:
  """
  touch annotation.bed12
  """
}
