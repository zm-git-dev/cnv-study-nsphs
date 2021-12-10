nextflow.enable.dsl=2

params.reference='/proj/sens2016007/nobackup/Reference'
params.bams="/proj/sens2016007/nobackup/bam-files/*.bam"
params.phenotypes = "/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData"
params.covariates="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv"
params.chromosomes=1..22

include {call_cnvs} from './cnv_calling/cnv_calling.nf'
include {create_matrix} from './cnv_matrix/cnv_matrix.nf'
include {gwas} from './gwas/gwas.nf'

workflow {
  call_cnvs(
    Channel.fromPath(params.bams),
    Channel.fromPath(params.reference)
  )
  create_matrix(
    call_cnvs.out.raw_variants,
    call_cnvs.out.qc_variants
  )
  gwas(
    create_matrix.out,
    Channel.fromPath(params.phenotypes),
    Channel.fromPath(params.covariates),
    Channel.of(params.chromosomes)
  )
}