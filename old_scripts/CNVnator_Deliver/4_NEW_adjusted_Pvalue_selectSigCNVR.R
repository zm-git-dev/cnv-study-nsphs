# During the summer project to lower the P-value threshold
#cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
# new p value: 438Proteins and 29400MainCNVRs
## Only keep P-value
rm(list=ls())
DataDir <- file.path("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables")
T_Pdf <- readRDS(file.path(DataDir, "T_All_P_only.rds"))

# for each biomarker, keep the loci with lowest P-value
Best_T_Pdf <- data.frame(P_value=apply(T_Pdf,1,min),Coor=colnames(T_Pdf)[apply(T_Pdf,1,which.min)])
#Format coordiante and biomarker name:
library(stringr)
Best_T_Pdf <- data.frame(Best_T_Pdf, Panel=str_split_fixed(rownames(Best_T_Pdf), "_", 3)[,1], BioMarker_name=str_split_fixed(rownames(Best_T_Pdf), "_", 3)[,3])
Coor <- Best_T_Pdf$Coor
Coor <- str_split_fixed(Coor, "[.]", 3)
Coor[,1] <- substring(Coor[,1], 2)
Coor <- apply(Coor, MARGIN = c(1,2), function(x) as.numeric(x) )
colnames(Coor) <- c("Chr", "Start", "End")
histData <- data.frame(Coor, Best_T_Pdf)
# Sort by Chr and Start:
histData_sort <- histData[order( histData$Chr, histData$Start), ]
histData_sort <- data.frame(number = c(1:nrow(histData_sort)),Best_T_Pdf)

#library("ggplot2")
# Manhattan plot for the highest values fo each biomarker: adapted from Nima Rafati nimarafati@gmail.com
# ggplot(data = histData_sort, aes(x = number, y = -log10(P_value), fill = factor(Panel))) +
#     geom_bar(stat = "identity") +
#     facet_grid(~Chr, scales = 'free_x', space = 'free_x', switch = 'x') +
#     theme_classic() +
#     theme(plot.title = element_text(size=22), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
#           strip.background = element_rect(colour = "white", fill = "grey"),panel.spacing = unit(0, "lines"), axis.text = element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12)) + scale_x_continuous(expand = c(0.1, 0.1)) +
#     xlab("Chromosome") +
#     ylim(0,60) +
#     ggtitle("Max -log10(P) for each biomarker") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     labs(fill="Panel") #+
#     #geom_text(data=histData_sort_filter,aes(x=number,y=-log10(P_value),label=BioMarker_name, size=5),vjust=0,hjust=0,angle=90,size=4)
# # Save best P for all biomarker manhattan plot as:
# aspect_ratio <- 2
# height <- 7
# ggsave("Manha_not_anno.png", height = height , width = height * aspect_ratio)
#setwd("/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS")
# keep regions pass adjusted P value for each of the 17 top biomarkers
#OLD: Pfilter <- 0.05/438/243987 # = 4.67e-10 0.05/438/29400=3.88e-9
AllValues <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/T_AllPvalues.rds")
Pfilter <- 0.05/438/29400 # 0.05/438/29400=3.88e-9
histData_sort_filter <- histData_sort[histData_sort$P_value <= 0.05/438/29400,] #dim 23X5
significantBio <- data.frame()
library(stringr)
WholeGenome_filter <- T_Pdf[rownames(histData_sort_filter),] #size: 23biomarkers 243,987AllPos with CNV
for(i in 1:nrow(WholeGenome_filter)) {
  OneBioMarker <- data.frame(P_value=WholeGenome_filter[i,][WholeGenome_filter[i,] <= Pfilter], Stat=AllValues[rownames(WholeGenome_filter)[i],][WholeGenome_filter[i,] <= Pfilter], Coor=colnames(WholeGenome_filter)[WholeGenome_filter[i,] <= Pfilter], Biomarker= rep(rownames(WholeGenome_filter)[i],sum(WholeGenome_filter[i,] <= Pfilter, na.rm = TRUE)))
  Coor <- OneBioMarker$Coor
  Coor <- str_split_fixed(Coor, "[.]", 3)
  Coor[,1] <- substring(Coor[,1], 2)
  Coor <- apply(Coor, MARGIN = c(1,2), function(x) as.numeric(x) )
  colnames(Coor) <- c("Chr", "Start", "End")
  OneBioMarker<- data.frame(Coor, OneBioMarker)
  significantBio <- rbind(significantBio, OneBioMarker)
}
write.csv(significantBio, file = "List_significantBiomarkers_withStat_changePvalue.csv")
#rm(AllValues)
# in Bash:
# change "," to tab, convert the csv file to bed format
sed 's/,/\t/g' List_significantBiomarkers_withStat_changePvalue.csv > List_significantBiomarkers_withStat_changePvalue.bed
# less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.bed
# remove header line:
sed -i '1d' List_significantBiomarkers_withStat_changePvalue.bed
# remove first column(coordiante): e.g. X1.161640580.161640780
cp List_significantBiomarkers_withStat_changePvalue.bed  tmp_List_significantBiomarkers_withStat_changePvalue.bed
#sed -i '1d' tmp_List_significantBiomarkers_withStat_changePvalue.bed
cut -f2- tmp_List_significantBiomarkers_withStat_changePvalue.bed > tmp_format_List_significantBiomarkers_withStat_changePvalue.bed
# number of significant associated biomarkers:
awk '{print $7}' tmp_format_List_significantBiomarkers_withStat_changePvalue.bed | sort | uniq | wc -l
# 23
module load bioinfo-tools BEDTools/2.27.1
bedtools sort -i tmp_format_List_significantBiomarkers_withStat_changePvalue.bed > SignificantCNPs.bed
awk '{OFS="\t"; print $1,$2,$3,$3-$2,$4,$7}' /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SignificantCNPs.bed > /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/List_significantBiomarkers_withStat_changePvalue_CNVnator_size.bed
bedtools merge -i SignificantCNPs.bed > SigCNVR_changePvalue.bed # this will be used to compare with Manta CNVRs
# add size:
awk '{OFS="\t"; print $0,$3-$2}' SigCNVR_changePvalue.bed > SigCNVR_size_changePvalue.bed # use this to intersect with two PacBio sequence available individual:
changeP_SigCNVR="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed" # compare with the two PacBio individual:

# NEW, use bedtools merge -o to optimize the merging and associated biomarker:
#sort -k1,1 -nk2,2 -nk3,3 tmp_format_List_significantBiomarkers_withStat_changePvalue.bed | bedtools merge -c 7,7 -o collapse,count | less -S
#bedtools merge -c 7,7 -o collapse,count -i tmp_format_List_significantBiomarkers_withStat_changePvalue.bed | less -S
#TODO: fix the collapse, now is just appending ALL######
#sort -k1,1 -nk2,2 -nk3,3 tmp_format_List_significantBiomarkers_withStat_changePvalue.bed | bedtools merge -c 7,7 -o collapse,count > SigCNVR_changePvalue_withBiomarker.bed
#cp SigCNVR_changePvalue_withBiomarker.bed  /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/.
#awk '{OFS="\t"; print $1,$2,$3,$3-$2,$4,$5}' SigCNVR_changePvalue_withBiomarker.bed | sort -k1,1 > resort_addSize_SigCNVR_changePvalue_withBiomarker.bed
#TODO: save to CSV and bed, generate a hg38 version also.
#cp resort_addSize_SigCNVR_changePvalue_withBiomarker.bed /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/.

# NEW fix collpase (e.g. uniq -c)
bedtools intersect -wao -a SigCNVR_size_changePvalue.bed -b tmp_format_List_significantBiomarkers_withStat_changePvalue.bed | bedtools groupby -g 1-4 -c 11 -o freqdesc > SigCNVR_size_changePvalue_withBiomarker.bed
cp SigCNVR_size_changePvalue_withBiomarker.bed /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/.
#convert to csv format in excel, add 分数点？

# Whole-genome Manhattan plot for the 23 top biomarkers
Manha23 <- data.frame()
library(stringr)
AllValues <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/T_AllPvalues.rds")
#WholeGenome_filter <- T_Pdf[rownames(histData_sort_filter),] #size: 17biomarkers 243987AllPos
# each of the 17 biomarkers, extract whole-genome wide CNVR results
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
#size: 23*243987
write.csv(Manha23, file = "Manha23.csv")
# Read Manha23.csv
#Manha23 <- read.csv(file="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/Manha23.csv", header=TRUE, sep=",")
colnames(Manha23) <- c("Coor", "CHR", "BP", "BP_end", "P", "Stat", "SNP", "Biomarker")
#colnames(Manha23) <- c("CHR", "BP", "BP_end", "P", "Stat", "SNP", "Biomarker")
sum(Manha23$P < 0.05/438/29400)
#458
#cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ManhaPlot
library(qqman)
png("Manha23_CoorHighlight.png", width=1080, height=720, res=120)
manhattan(Manha23, genomewideline = FALSE, suggestiveline = FALSE, highlight = Manha23$SNP[Manha23$P < 0.05/438/29400] )
dev.off()

#24Oct new: QQ-plot
png("Manha23_QQplot.png", width=1080, height=720, res=120)
qq(Manha23$P)
dev.off()
