process SAMTOOLS_SORT {
  label 'samtools'
  publishDir "${params.outdir}/star/${id}", mode: 'copy', pattern: "Aligned.sortedByCoord.out.bam"

  input:
  tuple val(id), path(bam), path(log_final), path(sj_tab)

  output:
  tuple val(id), path("Aligned.sortedByCoord.out.bam"), emit: sorted

  script:
  def mem_per_thread = (task.memory.toMega() / task.cpus).toInteger()
  """
  samtools sort -m ${mem_per_thread}M -@ ${task.cpus} \\
    -o Aligned.sortedByCoord.out.bam $bam > samtools_sort.log 2>&1
  """

  stub:
  """
  touch Aligned.sortedByCoord.out.bam
  """
}
