#!Rscript
library(foreach)
library(data.table)
library(doParallel)
library(magrittr)

Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1) %>%
  as.integer() %>%
  registerDoParallel()

load("!{phenotypes}")
AgeSex <- read.csv(file="!{covariates}", header=TRUE, sep=",")
CNVbed <- fread("!{cnv_matrix}")

LR_function <- function(BioMarker, CNVbed, AgeSex, Chrom) {
  #Select by chromosome:
  CNVbed <- CNVbed[chr==Chrom]
  cnv_coordinates <- CNVbed[, .(chr, start, end)]
  # Format data:
  CNVbed <- CNVbed[ , names(CNVbed) %in% rownames(BioMarker), with=FALSE]
  rownames(AgeSex) <- AgeSex$id
  AgeSex <- AgeSex[rownames(AgeSex) %in% colnames(CNVbed), ]
  BioMarker <- BioMarker[rownames(BioMarker) %in% colnames(CNVbed), ]
  CNVbed <- CNVbed[rowSums(is.na(CNVbed))/length(CNVbed)<0.1, ]
  BioMarker <- BioMarker[ ,colSums(is.na(BioMarker))/nrow(BioMarker)<0.1]
  CNVbed <- CNVbed[ ,order(names(CNVbed)), with=FALSE]
  CNVbed <- CNVbed[, !duplicated(names(CNVbed)), with=FALSE]
  AgeSex <- AgeSex[order(rownames(AgeSex)), ]
  BioMarker <- BioMarker[order(rownames(BioMarker)), ]
  Sex <- AgeSex[ ,2] # Population Sex
  Age <- AgeSex[ ,3] # Population Age
  
  foreach (i = 1:nrow(CNVbed), .combine=rbind) %dopar% {
    CNV <- t(CNVbed[i]) # index CNV
    foreach (i1 = 1:ncol(BioMarker), .combine=rbind) %dopar% {
      Y <- BioMarker[ ,i1] #index biomarkder
      model1 <- glm( Y ~ CNV + Sex + Age, family = gaussian, na.action = na.omit)
      coefs <- summary(model1)$coef[2, 1:4]

      cnv_coordinates[i, .(
        chr,
        start,
        end,
        marker = colnames(BioMarker)[i1],
        N = nobs(model1),
        beta = coefs[1],
        se = coefs[2],
        t = coefs[3],
        p = coefs[4]
      )]
    }
  }

}

gwas_results <- LR_function(pea_3, CNVbed, AgeSex, "!{chromosome}")
write.table(gwas_results, "!{chromosome}.glm", sep="\t", quote=F, row.names=F)