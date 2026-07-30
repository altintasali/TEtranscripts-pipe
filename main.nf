#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GUNZIP              } from './modules/gunzip'
include { DETECT_READ_LENGTH  } from './modules/detect_read_length'
include { STAR_INDEX          } from './modules/star_index'
include { GTF_TO_GENEPRED     } from './modules/gtf_to_genepred'
include { GENEPRED_TO_BED12   } from './modules/genepred_to_bed12'
include { STAR_ALIGN          } from './modules/star_align'
include { SAMTOOLS_SORT       } from './modules/samtools_sort'
include { SAMTOOLS_INDEX      } from './modules/samtools_index'
include { INFER_STRANDEDNESS  } from './modules/infer_strandedness'
include { DETERMINE_STRANDEDNESS } from './modules/determine_strandedness'
include { TECOUNT             } from './modules/tecount'
include { TETRANSCRIPTS       } from './modules/tetranscripts'
include { MULTIQC             } from './modules/multiqc'

workflow {
  checkParams()

  fasta_file   = file(params.fasta,   checkIfExists: true)
  gtf_file     = file(params.gtf,     checkIfExists: true)
  te_gtf_file  = file(params.te_gtf,  checkIfExists: true)
  py_script    = file("${projectDir}/scripts/determine_strandedness.py")

  // Parse samplesheet
  Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv(header: true, sep: ',')
    | map { parseSampleRow(it) }
    | set { samples_ch }

  // Sample metadata channel
  samples_ch
    | map { id, fq1, fq2, stranded_mode, condition ->
      [id, condition, stranded_mode]
    }
    | set { sample_meta_ch }

  // Gunzip references if needed
  GUNZIP(fasta_file, gtf_file, te_gtf_file)

  // Detect sjdbOverhang
  first_fastq = samples_ch
    | map { id, fq1, fq2, stranded_mode, condition -> fq1 }
    | first()

  DETECT_READ_LENGTH(first_fastq)

  sjdb_overhang = DETECT_READ_LENGTH.out
    | map { read_len ->
      params.sjdb_overhang == 'auto'
        ? read_len.toInteger() - 1
        : params.sjdb_overhang.toInteger()
    }

  // Build or use pre-built STAR index
  if (params.star_index) {
    star_idx_ch = Channel.fromPath(params.star_index, checkIfExists: true)
  } else {
    STAR_INDEX(GUNZIP.out.fasta, GUNZIP.out.gtf, sjdb_overhang)
    star_idx_ch = STAR_INDEX.out.index
  }

  // Annotation conversion for RSeQC
  GTF_TO_GENEPRED(GUNZIP.out.gtf)
  GENEPRED_TO_BED12(GTF_TO_GENEPRED.out.genepred)
  bed12_ch = GENEPRED_TO_BED12.out.bed12

  // Alignment
  align_input = samples_ch
    | map { id, fq1, fq2, stranded_mode, condition ->
      [id, fq1, fq2, stranded_mode, condition]
    }
    | combine(star_idx_ch.map { idx -> [idx] })

  STAR_ALIGN(align_input)

  // Sort and index BAMs
  SAMTOOLS_SORT(STAR_ALIGN.out.aligned)

  SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted)

  // Strandedness detection
  infer_input = SAMTOOLS_INDEX.out.indexed
    | join(sample_meta_ch)
    | combine(bed12_ch)

  INFER_STRANDEDNESS(infer_input)
  infer_ch = INFER_STRANDEDNESS.out.inferred

  // Split auto (needs determination) vs known
  infer_ch.branch {
    auto:  it[3] == 'auto'
    known: it[3] != 'auto'
  }.set { infer_split }

  // Determine strandedness for auto samples
  DETERMINE_STRANDEDNESS(infer_split.auto, py_script)

  det_stranded = DETERMINE_STRANDEDNESS.out.stranded
    | map { id, bam, stranded_file, condition ->
      [id, bam, stranded_file.text.trim(), condition]
    }

  // Use explicit value for known samples
  known_stranded = infer_split.known
    | map { id, bam, infer_txt, stranded_mode, condition ->
      [id, bam, stranded_mode, condition]
    }

  stranded_ch = det_stranded.mix(known_stranded)

  // TEcount
  tecount_input = stranded_ch
    | combine(GUNZIP.out.gtf)
    | combine(GUNZIP.out.te_gtf)

  TECOUNT(tecount_input)

  // TEtranscripts differential expression
  stranded_ch
    | map { id, bam, stranded_val, condition ->
      [condition, id, bam, stranded_val]
    }
    | groupTuple(by: 0)
    | toList()
    | flatMap { groups ->
      def contrasts = []
      for (int i = 0; i < groups.size(); i++) {
        for (int j = i + 1; j < groups.size(); j++) {
          def (c1_cond, c1_ids, c1_bams, c1_strands) = groups[i]
          def (c2_cond, c2_ids, c2_bams, c2_strands) = groups[j]

          def (treat_cond, treat_bams, treat_strands) =
            c1_cond < c2_cond
              ? [c2_cond, c2_bams, c2_strands]
              : [c1_cond, c1_bams, c1_strands]
          def (ctrl_cond, ctrl_bams, ctrl_strands) =
            c1_cond < c2_cond
              ? [c1_cond, c1_bams, c1_strands]
              : [c2_cond, c2_bams, c2_strands]

          def contrast_id = "${treat_cond}_vs_${ctrl_cond}"
          def all_strands = (treat_strands + ctrl_strands).unique()
          if (all_strands.size() > 1) {
            error "Inconsistent strandedness in contrast ${contrast_id}: ${all_strands}. " +
              "All samples in a contrast must have the same strandedness value."
          }
          contrasts << [contrast_id, treat_bams, ctrl_bams, all_strands[0]]
        }
      }
      return contrasts
    }
    | combine(GUNZIP.out.gtf)
    | combine(GUNZIP.out.te_gtf)
    | set { contrast_input_ch }

  TETRANSCRIPTS(contrast_input_ch)

  // MultiQC
  star_logs = STAR_ALIGN.out.aligned
    | map { id, bam, log_final -> log_final }
    | collect()

  rseqc_inputs = INFER_STRANDEDNESS.out.inferred
    | map { id, bam, infer_txt, mode, condition -> infer_txt }
    | collect()

  MULTIQC(star_logs, rseqc_inputs)
}

def checkParams() {
  def required = [
    samplesheet: params.samplesheet,
    fasta:       params.fasta,
    gtf:         params.gtf,
    te_gtf:      params.te_gtf
  ]
  def missing = required.findAll { k, v -> !v }.keySet()
  if (missing) {
    error "Missing required parameters: ${missing.join(', ')}"
  }
}

def parseSampleRow(row) {
  def id = row.sample?.trim()
  if (!id) error "Sample sheet row missing 'sample': ${row}"

  def fq1_str = row.fastq_1?.trim()
  if (!fq1_str) error "Sample '${id}' missing fastq_1"
  def fq1 = file(fq1_str, checkIfExists: true)

  def fq2_str = row.fastq_2?.trim()
  def fq2 = fq2_str ? file(fq2_str, checkIfExists: true) : []

  def stranded_mode = row.strandedness?.trim() ?: 'auto'
  stranded_mode = stranded_mode.toLowerCase()
  if (stranded_mode == 'unstranded') stranded_mode = 'no'

  def condition = row.condition?.trim() ?: ''

  [id, fq1, fq2, stranded_mode, condition]
}
