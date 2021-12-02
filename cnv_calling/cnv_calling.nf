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

process calculate_bins {
  cpus 1
  time '20m'
  input:
    path root
    path reference
  output:
    stdout
  shell:
    template 'calculate_bins.sh'
}

process partition {
  cpus 2
  time '40m'
  input:
    path root
    val bin_size
  output:
    path "${root}"
  shell:
    template 'partition.sh'
}

process call_variants {
  cpus 1
  time '1h'
  publishDir "cnv_calls", mode: 'symlink'
  input:
    path root
    val bin_size
  output:
    path '*_variants.txt'
  shell:
    template 'calling.sh'
}

process quality_control {
  cpus 1
  time '30m'
  publishDir "cnv_calls", mode: 'symlink'
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
  publishDir "cnv_calls", mode: 'symlink'
  stageInMode 'copy'
  
  input:
    path root
    path reference
  output:
    path '*_variants.txt'
  shell:
    template 'cnvnator.sh'
}

workflow {
  bams = Channel.fromPath(params.bams)
  extract_reads(bams)
  cnvnator(extract_reads.out, params.reference)
  quality_control(cnvnator.out)
}
