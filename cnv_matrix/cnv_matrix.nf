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

workflow create_matrix {
  take:
    raw_variants
    qc_variants
  main:
    filter_cnvs(raw_variants)
    make_windows(filter_cnvs.out.collect())
}