#cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SelectPacBio
#Population raw output"/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_population_bin10_rowname.rds"
# Sig CNVRs: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SigCNVR_size.bed

# module load R/3.4.3
# module load R_packages/3.4.3
#R:
setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SelectPacBio")
rm(list = ls())
CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_population_bin10_rowname.rds")
Sigbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SigCNVR_size.bed", sep = "\t")
colnames(Sigbed) <- c("CNVR_CHROM", "CNVR_Start", "CNVR_End", "CNVR_Size")
# for population matrix: handel NA and string type
# column names changed!
CNVbed_num <- data.frame(lapply(CNVbed, function(x) {as.numeric(as.character(x))}), check.names=FALSE)
#CNVbed_num <- data.frame(apply(CNVbed,2,function(x) {as.numeric(as.character(x))}), check.name=FALSE)
# fileter the Sig CNVRs by size:
#Sigbed <- Sigbed[Sigbed$CNVR_Size >= 2000, ]

SelectBed_df <- data.frame()
for(i in 1:nrow(Sigbed)) {
  print(paste0("Formatting CNVR NO.", i))
  dat <- subset(CNVbed_num, Chr==Sigbed[i,"CNVR_CHROM"] & Start>=Sigbed[i,"CNVR_Start"] & End<=Sigbed[i,"CNVR_End"])
  # calculate the mean copy number for each CNVR, individuals with one or more 'NA' will have NA copy number genotype in the CNVR
  dat_average <- apply(dat[,4:ncol(dat)], 2, function(x){mean(x)} ) #sum(is.na(SelectBed_df)):
  SelectBed <- data.frame(Sigbed[i, ], data.frame(t(dat_average), check.names=FALSE), check.names=FALSE)
  SelectBed_df <- rbind(SelectBed_df, SelectBed)
}
dim(SelectBed_df) # 26 1025
Ref_NA_count <- sort(apply(SelectBed_df[,5:ncol(SelectBed_df)], 2, function(x) {sum(x==2|is.na(x))} ))
write.table(Ref_NA_count, "Ref_NA_countt.txt", sep="\t")

Priority_df <- cbind(SelectBed_df[,1:4], SelectBed_df[,names(Ref_NA_count)])
write.table(Priority_df, "PriorityOrder_list_Manta.txt", sep="\t", row.names = F)
saveRDS(Priority_df, "PriorityOrder_list_Manta.rds")
