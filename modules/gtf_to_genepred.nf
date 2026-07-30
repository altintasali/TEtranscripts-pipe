process GTF_TO_GENEPRED {
  label 'genepred'

  input:
  path gtf

  output:
  path("annotation.genePred"), emit: genepred

  script:
  """
  gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \\
    $gtf annotation.genePred > gtf_to_genepred.log 2>&1
  """

  stub:
  """
  touch annotation.genePred
  """
}
