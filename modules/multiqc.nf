process MULTIQC {
  label 'multiqc'
  publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'multiqc_report.*'

  input:
  path star_logs
  path rseqc_results

  output:
  path("multiqc_report.html"), emit: html
  path("multiqc_report_data"),  emit: data

  script:
  """
  multiqc --force \\
    -o . -n multiqc_report \\
    $star_logs $rseqc_results \\
    > multiqc.log 2>&1
  """

  stub:
  """
  mkdir -p multiqc_report_data
  touch multiqc_report.html
  touch multiqc_report_data/multiqc_data.json
  """
}
