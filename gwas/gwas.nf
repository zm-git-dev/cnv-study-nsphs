nextflow.enable.dsl=2

process gwas_single_chromosome {
  cpus 4
  time '10h'
  beforeScript 'ml R_packages'
  publishDir "cnv_calls/gwas", mode: 'symlink'

  input:
    path cnv_matrix
    path phenotypes
    path covariates
    val chromosome
  output:
    path '*.glm'
  shell:
    template 'gwas.R'
}

workflow gwas { 
  take:
    cnv_matrix
    pheno
    covariates
    chromosomes
  main:
    gwas_single_chromosome(cnv_matrix, pheno, covariates, chromosomes)
  emit:
    gwas_single_chromosome.out
}

params.cnv_matrix = "/proj/sens2016007/nobackup/disentanglement/cnv_calls/matrix/cnv_matrix_collapsed.txt"
params.phenotypes = "/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData"
params.covariates="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv"
params.chromosomes=1..22

workflow {
  gwas(
    Channel.fromPath(params.cnv_matrix),
    Channel.fromPath(params.phenotypes),
    Channel.fromPath(params.covariates),
    Channel.of(params.chromosomes)
  )
}