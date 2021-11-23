# Scripts to include more biomerker (remove the  10% NA biomarker selection)
# 30 Oct to generate documentation for project conclusion
# module load R/3.4.3
# module load R_packages/3.4.3


library(data.table)
library(dplyr)

setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder_moreBiomarker")
Pdf <- data.frame()
for(i in 1:10) {
  fil <- paste(i, ".rds", sep = "")
  print(paste0("Reading chromosome:", i))
  temp <- readRDS(fil)
  Pdf <- rbind(Pdf, temp)
}
fil <- "AllPvalues1_10.rds"
saveRDS(Pdf,fil)

Pdf <- data.frame()
for(i in 11:22) {
  fil <- paste(i, ".rds", sep = "")
  print(paste0("Reading chromosome:", i))
  temp <- readRDS(fil)
  Pdf <- rbind(Pdf, temp)
}
fil <- "AllPvalues11_22.rds"
saveRDS(Pdf,fil)

Pdf <- data.frame()
temp <- readRDS("AllPvalues1_10.rds")
Pdf <- rbind(Pdf, temp)
temp <- readRDS("AllPvalues11_22.rds")
Pdf <- rbind(Pdf, temp)
# > dim(Pdf)
# [1] 108745    438
fil <- "AllPvalues.rds" #"
saveRDS(Pdf,fil)
T_Pdf <- t(Pdf)
fil <- "T_AllPvalues.rds"
saveRDS(T_Pdf,fil)

# select p value:
rm(list=ls())
T_Pdf <- readRDS("T_AllPvalues.rds")
#GET the last value of: Estimate/Std.Error/t-value/Pr(>|t|)
T_Pdf <- apply(T_Pdf, MARGIN = c(1,2), function(x) as.numeric(gsub("^.*/", "", x)) )
fil <- "T_All_P_only.rds"
saveRDS(T_Pdf,fil)
T_Pdf <- readRDS("T_All_P_only.rds")
# dim(T_Pdf) :438 108745
# number of regions? need for adjusted P-value threshold
#NOTE: wc -l /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/All.tmp.sort.merge.bed #=> 30598
#TODO: fix the NA for indexing and searching for lowest value
#Best_T_Pdf <- data.frame(P_value=apply(T_Pdf,1,min),Coor=colnames(T_Pdf)[apply(T_Pdf,1,which.min)])
#Best_T_Pdf <- data.frame(P_value=apply(T_Pdf,1,function(x){min(x, na.rm = TRUE)}),Coor=colnames(T_Pdf)[apply(T_Pdf,1,function(x){which (x == min(x, na.rm = TRUE))})])
Best_T_Pdf <- data.frame(P_value=apply(T_Pdf,1,min,na.rm=TRUE),Coor=colnames(T_Pdf)[apply(T_Pdf,1,which.min)])
library(stringr)
Best_T_Pdf <- data.frame(Best_T_Pdf, Panel=str_split_fixed(rownames(Best_T_Pdf), "_", 3)[,1], BioMarker_name=str_split_fixed(rownames(Best_T_Pdf), "_", 3)[,3], BioMarker_fullname=rownames(Best_T_Pdf))
Coor <- Best_T_Pdf$Coor
Coor <- str_split_fixed(Coor, "[.]", 3)
Coor[,1] <- substring(Coor[,1], 2)
Coor <- apply(Coor, MARGIN = c(1,2), function(x) as.numeric(x) )
colnames(Coor) <- c("Chr", "Start", "End")
histData <- data.frame(Coor, Best_T_Pdf)
# Sort by Chr and Start:
histData_sort <- histData[order( histData$Chr, histData$Start), ]
histData_sort <- data.frame(number = c(1:nrow(histData_sort)),Best_T_Pdf)
library(ggplot2)
# Manhattan plot for the highest values fo each biomarker: adapted from Nima Rafati nimarafati@gmail.com
ggplot(data = histData_sort, aes(x = number, y = -log10(P_value), fill = factor(Panel))) +
    geom_bar(stat = "identity") +
    facet_grid(~Chr, scales = 'free_x', space = 'free_x', switch = 'x') +
    theme_classic() +
    theme(plot.title = element_text(size=22), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          strip.background = element_rect(colour = "white", fill = "grey"),panel.spacing = unit(0, "lines"), axis.text = element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12)) + scale_x_continuous(expand = c(0.1, 0.1)) +
    xlab("Chromosome") +
    ylim(0,60) +
    ggtitle("Max -log10(P) for each biomarker") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Panel")

aspect_ratio <- 2
height <- 7
ggsave("Manha_not_anno_Manta.png", height = height , width = height * aspect_ratio)
#
AllValues <- readRDS("T_AllPvalues.rds")
Pfilter <- 0.05/438/29205 #0.05/397/20598 #update 30Oct
histData_sort_filter <- histData_sort[histData_sort$P_value <= 0.05/438/29205,] #dim 23*6 Old 0.05/397/30598; 25*9: 0.05/397/20598
significantBio <- data.frame()
library(stringr)
WholeGenome_filter <- T_Pdf[rownames(histData_sort_filter),] #[1]     23biomarkers 108745Pos
for(i in 1:nrow(WholeGenome_filter)) {
  OneBioMarker <- data.frame(P_value=WholeGenome_filter[i,][WholeGenome_filter[i,] <= Pfilter & !is.na(WholeGenome_filter[i,])], Stat=AllValues[rownames(WholeGenome_filter)[i],][WholeGenome_filter[i,] <= Pfilter & !is.na(WholeGenome_filter[i,])], Coor=colnames(WholeGenome_filter)[WholeGenome_filter[i,] <= Pfilter & !is.na(WholeGenome_filter[i,])], Biomarker= rep(rownames(WholeGenome_filter)[i],sum(WholeGenome_filter[i,] <= Pfilter, na.rm = TRUE)))
  Coor <- OneBioMarker$Coor
  Coor <- str_split_fixed(Coor, "[.]", 3)
  Coor[,1] <- substring(Coor[,1], 2)
  Coor <- apply(Coor, MARGIN = c(1,2), function(x) as.numeric(x) )
  colnames(Coor) <- c("Chr", "Start", "End")
  OneBioMarker<- data.frame(Coor, OneBioMarker)
  significantBio <- rbind(significantBio, OneBioMarker)
}
write.csv(significantBio, file = "List_significantBiomarkers_withStat_Manta.csv")
fil <- "List_significantBiomarkers_withStat_Manta.rds"
saveRDS(significantBio,fil)
write.table(significantBio, "List_significantBiomarkers_withStat_Manta.txt", sep="\t", col.names = F, row.names = F)

# in bash:
module load bioinfo-tools BEDTools/2.27.1
awk '{OFS="\t"; print $1,$2,$3,$3-$2,$4,$7}' List_significantBiomarkers_withStat_Manta.txt > List_significantBiomarkers_withStat_Manta_size.bed
# the SigCNVR for selecting individual:
bedtools sort -i List_significantBiomarkers_withStat_Manta_size.bed | bedtools merge -i stdin | awk '{OFS="\t";print $0,$3-$2}' > SigCNVR_size.bed
bedtools intersect -wao -a SigCNVR_size.bed -b List_significantBiomarkers_withStat_Manta_size.bed | bedtools groupby -g 1-4 -c 10 -o freqdesc > SigCNVR_size_withBiomarker_Manta.bed
#Manhattan plot:
# in R:

library(stringr)
AllValues <- readRDS("T_AllPvalues.rds")
#WholeGenome_filter <- T_Pdf[rownames(histData_sort_filter),] #size: 23biomarkers 108745Pos
# each of the 23 biomarkers, extract whole-genome wide CNVR results
Manha23 <- data.frame()
for(i in 1:nrow(WholeGenome_filter)) {
  OneBioMarker <- data.frame(P_value=WholeGenome_filter[i,], Stat=AllValues[rownames(WholeGenome_filter)[i],], Coor=colnames(WholeGenome_filter), Biomarker= rep(rownames(WholeGenome_filter)[i],length(WholeGenome_filter[i,]) ))
  Coor <- OneBioMarker$Coor
  Coor <- str_split_fixed(Coor, "[.]", 3)
  Coor[,1] <- substring(Coor[,1], 2)
  Coor <- apply(Coor, MARGIN = c(1,2), function(x) as.numeric(x) )
  colnames(Coor) <- c("Chr", "Start", "End")
  OneBioMarker<- data.frame(Coor, OneBioMarker)
  Manha23 <- rbind(Manha23, OneBioMarker)
}

# make plot:
#colnames(Manha23) <- c("Coor", "CHR", "BP", "BP_end", "P", "Stat", "SNP", "Biomarker")
fil <- "Manha23_Manta.rds"
saveRDS(Manha23,fil)
colnames(Manha23) <- c("CHR", "BP", "BP_end", "P", "Stat", "SNP", "Biomarker")
sum(Manha23$P < 0.05/438/29205, na.rm = TRUE) #[1] 130
#remove rows with NA P values:
Manha23_rmNA <- Manha23[!is.na(Manha23[,"P"]),]
dim(Manha23)
#[1] 2501135       7
dim(Manha23_rmNA)
#[1] 2328112       7

library(qqman)
png("Manha23_manta.png", width=1080, height=720, res=120)
manhattan(Manha23_rmNA, genomewideline = FALSE, suggestiveline = FALSE, highlight = Manha23_rmNA$SNP[Manha23_rmNA$P < 0.05/397/29205] )
dev.off()

#25Oct, QQplot
#Manha23 <- readRDS(Manha23_Manta.rds)
png("QQplot23_manta.png", width=1080, height=720, res=120)
qq(Manha23_rmNA$P)
dev.off()
