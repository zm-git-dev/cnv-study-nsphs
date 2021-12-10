#!Rscript
library(foreach)
library(data.table)

load("!{phenotypes}")
AgeSex <- read.csv(file="!{covariates}", header=TRUE, sep=",")
CNVbed <- fread("!{cnv_matrix}")

LR_function <- function(BioMarker, CNVbed, AgeSex, Chrom) {
  #Select by chromosome:
  CNVbed <- CNVbed[chr==Chrom]
  # Format data:
  CNVbed <- CNVbed[ , names(CNVbed) %in% rownames(BioMarker)]
  rownames(AgeSex) <- AgeSex$id
  AgeSex <- AgeSex[rownames(AgeSex) %in% colnames(CNVbed), ]
  BioMarker <- BioMarker[rownames(BioMarker) %in% colnames(CNVbed), ]
  CNVbed <- CNVbed[rowSums(is.na(CNVbed))/length(CNVbed)<0.1, ]
  BioMarker <- BioMarker[ ,colSums(is.na(BioMarker))/nrow(BioMarker)<0.1]
  CNVbed <- CNVbed[ ,order(names(CNVbed))]
  AgeSex <- AgeSex[order(rownames(AgeSex)), ]
  BioMarker <- BioMarker[order(rownames(BioMarker)), ]
  Sex <- AgeSex[ ,2] # Population Sex
  Age <- AgeSex[ ,3] # Population Age
  
  foreach (i = 1:nrow(CNVbed), .combine=rbind) %do% {
    CNV <- t(CNVbed[i, -(1:3)]) # index CNV
    foreach (i1 = 1:ncol(BioMarker), .combine=rbind) %do% {
      Y <- BioMarker[ ,i1] #index biomarkder
      model1 <- glm( Y ~ CNV + Sex + Age, family = gaussian, na.action = na.omit)
      coefs <- summary(model1)$coef[2, 1:4]

      CNVbed[i, .(chr, start, end, beta = coefs[1], se = coefs[2], t = coefs[3], p=coefs[4])]
    }
  }

}

LR_function(pea_3, CNVbed, AgeSex, "!{chromosome}")