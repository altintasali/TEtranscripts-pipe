process TETRANSCRIPTS {
  label 'tetranscripts'
  publishDir "${params.outdir}/tetranscripts/${contrast_id}", mode: 'copy'

  input:
  tuple val(contrast_id), path(treatment_bams), path(control_bams), val(stranded), path(gtf), path(te_gtf)

  output:
  tuple val(contrast_id), path("${contrast_id}.cntTable"), path("${contrast_id}_DESeq2.R"),
         path("${contrast_id}_gene_TE_analysis.txt"), path("${contrast_id}_sigdiff_gene_TE.txt"),
         emit: diffexp

  script:
  def treat_str = treatment_bams.collect { it.toString() }.join(' ')
  def ctrl_str = control_bams.collect { it.toString() }.join(' ')
  """
  TEtranscripts --format BAM --mode ${params.tetranscripts_mode} \\
    -t $treat_str \\
    -c $ctrl_str \\
    --GTF $gtf --TE $te_gtf \\
    --stranded $stranded \\
    --project $contrast_id \\
    --padj ${params.tetranscripts_padj} --foldchange ${params.tetranscripts_foldchange} \\
    --minread ${params.tetranscripts_minread} \\
    --outdir . \\
    ${params.tetranscripts_extra} \\
    > tetranscripts.log 2>&1
  """

  stub:
  """
  touch ${contrast_id}.cntTable
  touch ${contrast_id}_DESeq2.R
  touch ${contrast_id}_gene_TE_analysis.txt
  touch ${contrast_id}_sigdiff_gene_TE.txt
  """
}
