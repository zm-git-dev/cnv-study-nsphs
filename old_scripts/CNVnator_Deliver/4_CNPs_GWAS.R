# Genome-wide CNVnator:CNP association analysis in R:
## cat /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_LinearRegression.R
library(progress)
setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/P_tables")
#rm(list=ls())
## read pea_3 data:
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
AgeSex <- read.csv(file="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv", header=TRUE, sep=",")
#CNVbed <- read.table("CNV_1021.collapse.bed", header = TRUE, check.names=FALSE, sep = "\t") #the "-" is read as ".", use check.names
CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNVbed.rds")
BioCoor <- read.table("/proj/sens2016007/nobackup/Phenotypes/proteininfo_GRCh37_hg19_and_GRCh38_hg38_GENCODE.bed", header = FALSE, check.names=FALSE, sep = "\t")
LR_function <- function(BioMarker, BioCoor, CNVbed, AgeSex, Chrom) {
  #Select by chromosome:
  CNVbed <- CNVbed[CNVbed[,1]==Chrom,]
  #For whole genome, not selecting the cis-biomarkers
  #BioMarker <- BioMarker[,colnames(BioMarker) %in% BioCoor[BioCoor[,1]==Chrom,9]]
  # Format data:
  rownames(CNVbed) <- paste0(CNVbed$CHROM, ":", CNVbed$Start, "-", CNVbed$End)
  CNVbed <- CNVbed[ , names(CNVbed) %in% rownames(BioMarker)]
  rownames(AgeSex) <- AgeSex$id
  AgeSex <- AgeSex[rownames(AgeSex) %in% colnames(CNVbed), ]
  BioMarker <- BioMarker[rownames(BioMarker) %in% colnames(CNVbed), ]
  CNVbed[CNVbed<0] <- NA
  CNVbed <- CNVbed[rowSums(is.na(CNVbed))/length(CNVbed)<0.1, ]
  BioMarker <- BioMarker[ ,colSums(is.na(BioMarker))/nrow(BioMarker)<0.1]
  CNVbed <- CNVbed[ ,order(names(CNVbed))]
  AgeSex <- AgeSex[order(rownames(AgeSex)), ]
  BioMarker <- BioMarker[order(rownames(BioMarker)), ]
  pb <- progress_bar$new(
  format = "  Processing :what [:bar] :percent eta: :eta",
  clear = FALSE, total = nrow(CNVbed), width = 60)
  pb$tick(0)
  #for(i in 1:ncol(BioMarker)) {
  for(i in 1:nrow(CNVbed)) {
    CNV <- t(CNVbed[i, ]) # index CNV
    #TODO: remove the output file in case of append.
    #colnames( sort.new_testPea_3 )[i] # Get column's name, save to each file:
    #Y <- BioMarker[ ,i] #index biomarkder
    Sex <- AgeSex[ ,2] # Population Sex
    Age <- AgeSex[ ,3] # Population Age

    # in the second loop, save them to a table and save the whole table/list? in the end of the loop?
    #for(i1 in 1:nrow(CNVbed)) {
    for(i1 in 1:ncol(BioMarker)) {
      #CNV <- t(CNVbed[i1, ]) # index CNV
      Y <- BioMarker[ ,i1] #index biomarkder
      model1 <- glm( Y ~ CNV + Sex + Age, family = gaussian, na.action = na.omit)
      #TODO: save it to a list, when this loop is finish, save this list as a column to the output table:
      #write(summary(model1)$coef[2,],file=paste0(colnames( sort.new_testPea_3 )[i], ".txt"),append=TRUE)
      if (i1 == 1) {
        result <- paste(c(summary(model1)$coef[2, 1:4]), collapse='/' )
        #result <- summary(model1)$coef[2, 4]
      } else {
        result <- c(result, paste(c(summary(model1)$coef[2, 1:4]), collapse='/' ))
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
    pb$tick(tokens = list(what = Chrom))
  }
  LRresult <- t(LRresult)
  colnames(LRresult) <- colnames(BioMarker)
  fil <- paste(Chrom, ".rds", sep = "")
  saveRDS(LRresult,fil)
  # TODO: save the result by name of chromosome:
  printLog <- paste("LR analyses for CNVRs in chromosome:", Chrom, "finished. Result saved as", fil, sep = " ")
  return(printLog)
}

# split the association analysis by chromosome:
cd /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun
for i in {1..22}
do
   cat > CNV_Pea3_Chr${i}.sh <<- EOM
#!/bin/bash -l
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 9:00:00
#SBATCH -J CNV_Pea3_Chr${i}

module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
R CMD BATCH --no-restore --no-save CNV_Pea3_Chr${i}.R
EOM
done

for i in {1..22}
do
   cat > CNV_Pea3_Chr${i}.R <<- EOM
source("/home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_LinearRegression.R")
LR_function(pea_3, BioCoor, CNVbed, AgeSex, ${i})
EOM
done

for each in *.sh; do echo "sbatch ${each}"; done

ls /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/
## After 1st rerun: need to rerun 2,9
for i in {2,9}
do
   cat > CNV_Pea3_Chr${i}.sh <<- EOM
#!/bin/bash -l
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 19:00:00
#SBATCH -J CNV_Pea3_Chr${i}

module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
R CMD BATCH --no-restore --no-save CNV_Pea3_Chr${i}.R
EOM
done

# Rbind the P-values table for 22 chromosomes:
touch /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun/P_rbind.sh
vi /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun/P_rbind.sh
###########
#!/bin/bash -l
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5:00:00
#SBATCH -J P_rbind

module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
R CMD BATCH --no-restore --no-save /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun/P_rbind.R

# R scripts
touch /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun/P_rbind.R
vi /home/zhiwei94/CNVscripts/CNV_GWASfolder/WholeG_splitRun/P_rbind.R
#########
setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables")
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


# read all 22 rds and combine as one table:
module load R/3.4.3
module load R_packages/3.4.3
setwd("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables")
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
fil <- "AllPvalues.rds" #"
saveRDS(Pdf,fil)
T_Pdf <- t(Pdf)
fil <- "T_AllPvalues.rds"
saveRDS(T_Pdf,fil)

rm(list=ls())
T_Pdf <- readRDS("T_AllPvalues.rds")
#GET the last value of: Estimate/Std.Error/t-value/Pr(>|t|)
T_Pdf <- apply(T_Pdf, MARGIN = c(1,2), function(x) as.numeric(gsub("^.*/", "", x)) )
fil <- "T_All_P_only.rds"
saveRDS(T_Pdf,fil)

#cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables
## Only keep P-value
rm(list=ls())
T_Pdf <- readRDS("T_AllPvalues.rds")
#GET the last value of: Estimate/Std.Error/t-value/Pr(>|t|)
T_Pdf <- apply(T_Pdf, MARGIN = c(1,2), function(x) as.numeric(gsub("^.*/", "", x)) )
fil <- "T_All_P_only.rds"
saveRDS(T_Pdf,fil)
# read file
T_Pdf <- readRDS("T_All_P_only.rds")
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
histData_sort <- data.frame(number = c(1:nrow(histData_sort)),histData_sort)
library("ggplot2")

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
    labs(fill="Panel") #+
    #geom_text(data=histData_sort_filter,aes(x=number,y=-log10(P_value),label=BioMarker_name, size=5),vjust=0,hjust=0,angle=90,size=4)
# Save best P for all biomarker manhattan plot as:
aspect_ratio <- 2
height <- 7
ggsave("Manha_not_anno.png", height = height , width = height * aspect_ratio)
#setwd("/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS")
# keep regions pass adjusted P value for each of the 17 top biomarkers
AllValues <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/T_AllPvalues.rds")
Pfilter <- 0.05/438/243987
histData_sort_filter <- histData_sort[histData_sort$P_value <= 0.05/438/243987,] #dim 17X8
significantBio <- data.frame()
library(stringr)
WholeGenome_filter <- T_Pdf[rownames(histData_sort_filter),] #size: 17biomarkers 243987AllPos
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
write.csv(significantBio, file = "List_significantBiomarkers_withStat.csv")
# in Bash:
# change "," to tab, convert the csv file to bed format
sed 's/,/\t/g' List_significantBiomarkers_withStat.csv > List_significantBiomarkers_withStat.bed
# less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.bed
# remove header line:
sed -i '1d' List_significantBiomarkers_withStat.bed
# remove first column(coordiante): e.g. X1.161640580.161640780
cp tmpList_significantBiomarkers_withStat.bed  tmp2List_significantBiomarkers_withStat.bed
sed -i '1d' tmp2List_significantBiomarkers_withStat.bed
cut -f2- tmp2List_significantBiomarkers_withStat.bed > Format_tmp2List_significantBiomarkers_withStat.bed
module load bioinfo-tools BEDTools/2.27.1
bedtools sort -i Format_tmp2List_significantBiomarkers_withStat.bed > SignificantCNPs.bed
bedtools merge -i SignificantCNPs.bed # this will be used to compare with Manta CNVRs

# Whole-genome Manhattan plot for the 17 top biomarkers
Manha17 <- data.frame()
library(stringr)
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
  Manha17 <- rbind(Manha17, OneBioMarker)
}
#size: 17*243987=4 147 779
write.csv(Manha17, file = "Manha17.csv")
# Read Manha17.csv
Manha17 <- read.csv(file="Manha17.csv", header=TRUE, sep=",")
colnames(Manha17) <- c("Coor", "CHR", "BP", "BP_end", "P", "Stat", "SNP", "Biomarker")
sum(Manha17$P < 0.05/438/243987)
#382
library(qqman)
png("Manha17_CoorHighlight.png", width=1080, height=720, res=120)
manhattan(Manha17, genomewideline = FALSE, suggestiveline = FALSE, highlight = Manha17$SNP[Manha17$P < 0.05/438/243987] )
dev.off()
