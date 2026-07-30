process STAR_ALIGN {
  label 'star_align'
  publishDir "${params.outdir}/star/${id}", mode: 'copy', pattern: 'star_out/**/*'

  input:
  tuple val(id), path(fq1), path(fq2), val(stranded_mode), val(condition), path(star_idx)

  output:
  tuple val(id), path("star_out/${id}/Aligned.out.bam"), path("star_out/${id}/Log.final.out"),
         path("star_out/${id}/SJ.out.tab"), emit: aligned

  script:
  def read_command = ''
  if (fq1.name.endsWith('.gz')) {
    read_command = '--readFilesCommand zcat'
  }

  def fq_str = fq1
  if (fq2) {
    fq_str = "$fq1 $fq2"
  }

  """
  mkdir -p star_out/${id}
  STAR --runThreadN ${task.cpus} \\
    --genomeDir $star_idx \\
    --readFilesIn $fq_str \\
    $read_command \\
    --outFileNamePrefix star_out/${id}/ \\
    --outSAMtype BAM Unsorted \\
    ${params.star_extra} \\
    > star_out/${id}/align.log 2>&1
  """

  stub:
  """
  mkdir -p star_out/${id}
  touch star_out/${id}/Aligned.out.bam
  touch star_out/${id}/Log.final.out
  touch star_out/${id}/SJ.out.tab
  """
}
