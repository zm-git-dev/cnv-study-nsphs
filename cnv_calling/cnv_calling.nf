nextflow.enable.dsl=2

// Parameters
params.reference='/proj/sens2016007/nobackup/Reference'
params.bams="/proj/sens2016007/nobackup/bam-files/*.bam"

process extract_reads {
  cpus 4
  time '2h'
  input:
    path bam
  output:
    path '*.root'
  shell:
    template 'extract_reads.sh'
}

process quality_control {
  cpus 1
  time '30m'
  publishDir "cnv_calls/qc", mode: 'copy'
  beforeScript 'ml R_packages'

  input:
    path variants
  output:
    path '*_variants_qc.bed'
  shell:
    template 'qc.R'
}

process cnvnator {
  cpus 4
  time '4h'
  publishDir "cnv_calls/raw", mode: 'copy'
  stageInMode 'copy'
  
  input:
    path root
    path reference
  output:
    path '*_variants.txt'
  shell:
    template 'cnvnator.sh'
}

workflow call_cnvs {
  take:
    bams
    reference
  main:
    extract_reads(bams)
    cnvnator(extract_reads.out, reference)
    quality_control(cnvnator.out)
  emit:
    raw_variants = cnvnator.out
    qc_variants = quality_control.out
}

workflow {
  bams = Channel.fromPath(params.bams)
  call_cnvs(bams, params.reference)
}