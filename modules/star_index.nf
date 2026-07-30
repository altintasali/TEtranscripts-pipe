process STAR_INDEX {
  label 'star_index'
  publishDir "${params.outdir}/star_index", mode: 'copy', pattern: 'star_index/**/*'

  input:
  path fasta
  path gtf
  val  sjdb_overhang

  output:
  path("star_index"), emit: index

  script:
  """
  mkdir -p star_index
  STAR --runMode genomeGenerate \\
    --genomeDir star_index \\
    --genomeFastaFiles $fasta \\
    --sjdbGTFfile $gtf \\
    --sjdbOverhang $sjdb_overhang \\
    --runThreadN ${task.cpus} \\
    ${params.star_index_extra} \\
    > star_index/index.log 2>&1
  """

  stub:
  """
  mkdir -p star_index
  touch star_index/SA star_index/SAindex star_index/Genome
  """
}
