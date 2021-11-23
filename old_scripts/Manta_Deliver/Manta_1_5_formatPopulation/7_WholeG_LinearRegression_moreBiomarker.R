#rm(list=ls())
setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder_moreBiomarker")
rm(list=ls())
## read pea_3 data:
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
AgeSex <- read.csv(file="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv", header=TRUE, sep=",")
CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_population_bin10_rowname.rds")
##BioCoor <- read.table("/proj/sens2016007/nobackup/Phenotypes/proteininfo_GRCh37_hg19_and_GRCh38_hg38_GENCODE.bed", header = FALSE, check.names=FALSE, sep = "\t")

LR_function <- function(BioMarker,CNVbed, AgeSex, Chrom) {
#LR_function <- function(BioMarker, BioCoor, CNVbed, AgeSex, Chrom) {
  #Select by chromosome:
  CNVbed <- CNVbed[CNVbed[,1]==Chrom,]
  #For whole genome, not selecting the cis-biomarkers
  #BioMarker <- BioMarker[,colnames(BioMarker) %in% BioCoor[BioCoor[,1]==Chrom,9]]
  # Format data:
  #rownames(CNVbed) <- paste0(CNVbed$Chrom, ":", CNVbed$Start, "-", CNVbed$End)
  CNVbed <- CNVbed[ , names(CNVbed) %in% rownames(BioMarker)] #BioMarker <- pea_3
  rownames(AgeSex) <- AgeSex$id
  AgeSex <- AgeSex[rownames(AgeSex) %in% colnames(CNVbed), ]
  BioMarker <- BioMarker[rownames(BioMarker) %in% colnames(CNVbed), ]
  #CNVbed[CNVbed<0] <- NA
  CNVbed <- CNVbed[rowSums(is.na(CNVbed))/length(CNVbed)<0.1, ] # over 10% high quality
  #CNVbed <- CNVbed[rowSums(CNVbed=="2") < length(CNVbed) - 1, ] # at least one sample not 2 copies(default, non-detected genotype), case with one not 2copies error
  #CNVbed <- CNVbed[rowSums(is.na(CNVbed)) < length(CNVbed) - 1, ]
  #BioMarker <- BioMarker[ ,colSums(is.na(BioMarker))/nrow(BioMarker)<0.1]
  CNVbed <- CNVbed[ ,order(names(CNVbed))]
  AgeSex <- AgeSex[order(rownames(AgeSex)), ]
  BioMarker <- BioMarker[order(rownames(BioMarker)), ]
  Sex <- AgeSex[ ,2] # Population Sex
  Age <- AgeSex[ ,3] # Population Age
  # pb <- progress_bar$new(
  # format = "  Processing :what [:bar] :percent eta: :eta",
  # clear = FALSE, total = nrow(CNVbed), width = 60)
  # pb$tick(0)
  # the above scripts are fine!
  #for(i in 1:ncol(BioMarker)) {
  for(i in 1:nrow(CNVbed)) {
    CNV <- t(CNVbed[i, ]) # index CNV ?
    #print(i)
    print(paste0("Analysing:",rownames(CNVbed)[i]))
    #print(summary(t(CNVbed[i, ])))
    #TODO: remove the output file in case of append.
    #colnames( sort.new_testPea_3 )[i] # Get column's name, save to each file:
    #Y <- BioMarker[ ,i] #index biomarkder
    #Sex <- AgeSex[ ,2] # Population Sex
    #Age <- AgeSex[ ,3] # Population Age

    # in the second loop, save them to a table and save the whole table/list? in the end of the loop?
    #for(i1 in 1:nrow(CNVbed)) {
    for(i1 in 1:ncol(BioMarker)) {
      #CNV <- t(CNVbed[i1, ]) # index CNV
      Y <- BioMarker[ ,i1] #index biomarkder
      #Y_name <- colnames(BioMarker)[i1]
      #glm_dt <- data.frame(Y,CNV,Sex,Age)
      #names(glm_dt)[names(glm_dt) == "Y"] <- Y_name
      #test_glm_dt <<- glm_dt
      #print(unique(glm_dt[is.na(glm_dt[,1]) == 0,2]))
      #if (length(unique(glm_dt[is.na(glm_dt[,1]) == 0,2])) == 1) {
      if (length(unique(na.omit(CNV[is.na(Y) == 0]))) == 1) {
        result_one <- NA
      } else {
        model1 <- glm( Y ~ CNV + Sex + Age, family = gaussian, na.action = na.omit) ##!!! bug # remove regions with all 2 copies genotype
        result_one <- paste(c(summary(model1)$coef[2, 1:4]), collapse='/' )
      }
      #TODO: save it to a list, when this loop is finish, save this list as a column to the output table:
      #write(summary(model1)$coef[2,],file=paste0(colnames( sort.new_testPea_3 )[i], ".txt"),append=TRUE)
      if (i1 == 1) {
        result <- result_one
        #result <- summary(model1)$coef[2, 4]
      } else {
        result <- c(result, result_one)
        #result <- c(result, summary(model1)$coef[2, 4])
      }
    }
    result.matrix <- as.matrix(result)
    colnames(result.matrix) <- rownames(CNVbed)[i]
    #rownames(result.matrix) <- rownames(CNVbed)
    if (i == 1) {
      #col.name <- colnames(BioMarker)[i]
      LRresult <- result.matrix
    } else {
      LRresult <- data.frame(LRresult, result.matrix)
    }
    #pb$tick(tokens = list(what = Chrom))
  }
  LRresult <- t(LRresult)
  colnames(LRresult) <- colnames(BioMarker)
  fil <- paste(Chrom, ".rds", sep = "")
  saveRDS(LRresult,fil)
  # TODO: save the result by name of chromosome:
  printLog <- paste("LR analyses for CNVRs in chromosome:", Chrom, "finished. Result saved as", fil, sep = " ")
  return(printLog)
}
