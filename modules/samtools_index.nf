process SAMTOOLS_INDEX {
  label 'samtools_index'
  publishDir "${params.outdir}/star/${id}", mode: 'copy', pattern: "Aligned.sortedByCoord.out.bam.bai"

  input:
  tuple val(id), path(sorted_bam)

  output:
  tuple val(id), path(sorted_bam), path("${sorted_bam}.bai"), emit: indexed

  script:
  """
  samtools index -@ ${task.cpus} $sorted_bam > samtools_index.log 2>&1
  """

  stub:
  """
  touch $sorted_bam.bai
  """
}
