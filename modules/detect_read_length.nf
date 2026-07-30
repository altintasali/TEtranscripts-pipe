process DETECT_READ_LENGTH {
  label 'minimal'

  input:
  path fastq

  output:
  stdout

  script:
  """
  zcat -f $fastq | awk 'NR%4==2{if(length>m)m=length}END{print m}'
  """

  stub:
  """
  echo "101"
  """
}
