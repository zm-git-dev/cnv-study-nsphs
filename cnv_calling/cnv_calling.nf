nextflow.enable.dsl=2

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

workflow {
  bams = Channel.fromPath("/proj/sens2016007/nobackup/bam-files/*.bam")
  extract_reads(bams)

}