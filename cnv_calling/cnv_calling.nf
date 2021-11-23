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
    path 'out'
  script:
    template 'extract_reads.sh'
}

process calculate_bins {
  cpus 1
  time '20:00'

  input:
    path bam, root
  output:
    
  script:
    template 'calculate_bins.sh'
}

workflow {
  bams = Channel.fromPath(params.bams)
  extract_reads(bams)

}