nextflow.enable.dsl=2

// Parameters
params.reference='/proj/sens2016007/nobackup/Reference'
params.bams="/proj/sens2016007/nobackup/bam-files/*.bam"

process extract_reads {
  cpus 2
  time '1h'
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

process call {
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
  time '10m'
  publishDir "cnv_calls", mode: 'symlink'
  conda 'r-data.table r-tidyverse'
  input:
    path variants
  output
    path '*_variants_qc.bed'
  shell:
    template 'qc.R'
}

workflow {
  bams = Channel.fromPath(params.bams)
  extract_reads(bams)
  calculate_bins(extract_reads.out, params.reference)
  partition(extract_reads.out, calculate_bins.out)
  call(partition.out, calculate_bins.out)
  quality_control(call.out)
}