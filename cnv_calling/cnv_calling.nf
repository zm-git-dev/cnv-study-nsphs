nextflow.enable.dsl=2

// Parameters
params.reference='/proj/sens2016007/nobackup/Reference'
params.bams="/proj/sens2016007/nobackup/bam-files/*.bam"

process extract_reads {
  cpus 2
  time '1:00:00'
  input:
    path bam
  output:
    path 'root'
  shell:
    template 'extract_reads.sh'
}

process calculate_bins {
  cpus 1
  time '20:00'
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
  time '40:00'
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
  time '1:00:00'
  input:
    path root
    val bin_size
  output:
    path 'variants.txt'
  shell:
    template 'calling.sh'
}

process quality_control {
  
}

workflow {
  bams = Channel.fromPath(params.bams)
  extract_reads(bams)
  calculate_bins(extract_reads.out, params.reference)
  calculate_bins.out.view()
}