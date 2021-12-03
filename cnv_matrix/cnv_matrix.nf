nextflow.enable.dsl=2

process filter_cnvs {
  cpus 1
  time '10m'
  beforeScript 'ml R_packages'

  input:
    path variants
  output:
    path "*_filtered.bed"
  shell:
    template 'filter_variants.R'
}

process make_windows {
  cpus 4
  time '30m'

  input:
    path variants
  output:
    path 'windows.bed'
  shell:
    template 'makewindows.sh'
}

process align_cnvs {
  cpus 2
  time '1h'

  input:
    path variants
    path windows
    val sample_size
  output:
    path '*_200bp.bed'
  shell:
    template 'align_cnvs.bed'
}


workflow create_matrix {
  take:
    raw_variants
    qc_variants
  main:
    filter_cnvs(raw_variants)
    filtered_cnvs = filter_cnvs.out.collect()
    num_samples = filtered_cnvs.size()
    make_windows(filtered_cnvs)
  emit:
    make_windows.out
}

params.raw_variants="/proj/sens2016007/nobackup/disentanglement/cnv_calling/cnv_calls/raw/*"
params.qc_variants="/proj/sens2016007/nobackup/disentanglement/cnv_calling/cnv_calls/qc/*"
workflow {
  create_matrix(params.raw_variants, params.qc_variants)
  create_matrix.out.view()
}