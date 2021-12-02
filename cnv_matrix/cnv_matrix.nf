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

workflow create_matrix {
  take:
    raw_variants
    qc_variants
  main:
    filter_cnvs(raw_variants)
}