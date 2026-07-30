process GUNZIP {
  label 'minimal'

  input:
  path fasta
  path gtf
  path te_gtf

  output:
  path("fasta/*"),  emit: fasta
  path("gtf/*"),    emit: gtf
  path("te_gtf/*"), emit: te_gtf

  script:
  def gunzip_if_needed = { path, subdir ->
    if (path.name.endsWith('.gz')) {
      """
      mkdir -p ${subdir}
      gunzip -c ${path} > ${subdir}/${path.name.substring(0, path.name.length() - 3)}
      """
    } else {
      """
      mkdir -p ${subdir}
      ln -sf ${path} ${subdir}/${path.name}
      """
    }
  }

  """
  ${gunzip_if_needed(fasta, 'fasta')}
  ${gunzip_if_needed(gtf, 'gtf')}
  ${gunzip_if_needed(te_gtf, 'te_gtf')}
  """

  stub:
  def stub_link = { path, subdir ->
    def out = path.name.endsWith('.gz') ? path.name.substring(0, path.name.length() - 3) : path.name
    """
    mkdir -p ${subdir}
    touch ${subdir}/${out}
    """
  }
  """
  ${stub_link(fasta, 'fasta')}
  ${stub_link(gtf, 'gtf')}
  ${stub_link(te_gtf, 'te_gtf')}
  """
}
