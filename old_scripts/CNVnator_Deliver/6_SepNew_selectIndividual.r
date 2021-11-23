#cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/SelectPacBio
#changeP_SigCNVR="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
# population raw output: less -S  /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNV_1021.collapse.bed

# module load R/3.4.3
# module load R_packages/3.4.3
#R:
setwd("/proj/sens2016007/nobackup/Zhiwei/SummerProject/SelectPacBio")
rm(list = ls())
CNVbed <- readRDS("/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNVbed.rds")
Sigbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed", sep = "\t")
colnames(Sigbed) <- c("CNVR_CHROM", "CNVR_Start", "CNVR_End", "CNVR_Size")
# for population matrix: set negative value(low quality) to NA:
CNVbed[CNVbed<0] <- NA
# fileter the Sig CNVRs by 2kb
Sigbed <- Sigbed[Sigbed$CNVR_Size >= 2000, ]

SelectBed_df <- data.frame()
for(i in 1:nrow(Sigbed)) {
  print(paste0("Formatting CNVR NO.", i))
  dat <- subset(CNVbed, CHROM==Sigbed[i,"CNVR_CHROM"] & Start>=Sigbed[i,"CNVR_Start"] & End<=Sigbed[i,"CNVR_End"])
  # calculate the mean copy number for each CNVR, individuals with one or more 'NA' will have NA copy number genotype in the CNVR
  dat_average <- apply(dat[,4:ncol(dat)], 2, function(x){mean(x)} ) #sum(is.na(SelectBed_df)): with na.omit=TRUE:145 with out: 925
  SelectBed <- data.frame(Sigbed[i, ], data.frame(t(dat_average), check.names=FALSE), check.names=FALSE)
  SelectBed_df <- rbind(SelectBed_df, SelectBed)
}
# exclude MHC regions:
# vector with sample names and number of Ref/NA regions, sort from low(best) to high value.
noMHC <- sort(apply(SelectBed_df[c(1:17,26),5:ncol(SelectBed_df)], 2, function(x) {sum(x==2|is.na(x))} ))
write.table(noMHC, "noMHC_Ref_NA_countt.txt", sep="\t")
# noMHC <- sort(apply(SelectBed_df[c(1:17,26),5:ncol(SelectBed_df)], 2, function(x) {sum(x==2|is.na(x))} ))[1:64]
# withMHC <- sort(apply(SelectBed_df[,5:ncol(SelectBed_df)], 2, function(x) {sum(x==2|is.na(x))} ))[1:64]
# check_both <- names(noMHC) %in% names(withMHC)
# names(noMHC)[check_both]
# prioritisation list:
#Priority_df <- cbind(SelectBed_df[,1:4], SelectBed_df[names(noMHC)[check_both]])
Priority_df <- cbind(SelectBed_df[,1:4], SelectBed_df[,names(noMHC)])
# NA_2_count <- c("NO","MHC","Ref/NA","Count",noMHC)
# names(NA_2_count)[1:4] <-  c("CNVR_CHROM", "CNVR_Start", "CNVR_End", "CNVR_Size")
# test_df <- rbind(NA_2_count, Priority_df)
# save the prioritisation list:
write.table(Priority_df, "PriorityOrder_list.txt", sep="\t", row.names = F)
saveRDS(Priority_df, "PriorityOrder_list.rds")
# reorder whole population:
#order_SelectBed_df <- cbind(SelectBed_df[,1:4], SelectBed_df[,order])
#
