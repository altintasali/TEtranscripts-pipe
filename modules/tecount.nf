process TECOUNT {
  label 'tecount'
  publishDir "${params.outdir}/tecount/${id}", mode: 'copy', pattern: "${id}.cntTable"

  input:
  tuple val(id), path(bam), val(stranded), val(condition), path(gtf), path(te_gtf)

  output:
  tuple val(id), path("${id}.cntTable"), emit: cnt_table

  script:
  """
  mkdir -p tecount_out
  TEcount --format BAM --mode ${params.tetranscripts_mode} \\
    -b $bam \\
    --GTF $gtf --TE $te_gtf \\
    --stranded $stranded \\
    --project $id \\
    --outdir tecount_out \\
    ${params.tetranscripts_extra} \\
    > tecount.log 2>&1
  cp tecount_out/${id}.cntTable .
  """

  stub:
  """
  echo -e 'gene\\t${id}' > ${id}.cntTable
  """
}
