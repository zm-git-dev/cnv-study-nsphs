 # Compare CNVnator significant CNVRs with Manta population CNVRs:
## CNVnator folder: /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/StatFolder
## Manta folder: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/
CNVnatorDUP=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/StatFolder/All.tmp.sort.DUP.merge.bed
CNVnatorDEL=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/StatFolder/All.tmp.sort.DEL.merge.bed
# Signifiant Positions: CNVnator GWAS:/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS
SigBED=/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged.bed
MantaDUP=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/All.tmp.sort.DUP.mergeCount.bed
MantaDEL=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/All.tmp.sort.DEL.mergeCount.bed
module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -f 0.8 -wao -a $SigBED -b $MantaDEL | less -S
19      41381725        41387525        19      41344495        41387879        523     43384   5800
19      51508940        51510740        19      51508588        51511299        10      2711    1800
19      54555500        54560700        19      54555339        54560900        306     5561    5200

# add 4th column: size of the CNVRs: /proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged.bed
awk -v OFS='\t' '{print $0, $3-$2}' /proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged.bed > /proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed
SigBED_size=/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed

bedtools intersect -f 0.8 -wao -a $SigBED_size -b All.filter.tmp.bed | less -S # tried record overlap?
bedtools intersect -f 0.8 -wao -a $SigBED_size -b All.filter.extra20Mb.tmp.bed | less -S
sort -nk9,9 CNVsig_Manta100M.bed | less -S
sort  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/StatFinalPre/CNVsig_Manta100M_sorted.bed | uniq -c | sort -nk2,2 -nk10,10  | less -S

#6June: filter large CNVs:
cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/StatFinalPre
awk -v OFS="\t" '($4 <= 100000000 ){print $1, $2, $3, $15, $4}' ../../*PASS.bed > All.filter.tmp.bed #filter by size
sed  -i -e '/MT/d' All.filter.tmp.bed
# 10June extra filter Manta by 20Mb
awk -v OFS="\t" '($5 <= 20000000){print $0}' All.filter.tmp.bed > All.filter.extra20Mb.tmp.bed
sort -k 1,1 -k2,2n All.filter.extra20Mb.tmp.bed > All.filter.extra20Mb.sort.bed
bedtools merge -c 1 -o count -i All.filter.extra20Mb.sort.bed > Manta20Mb.sort.merge.bed
awk -v OFS='\t' '{print $0, $3-$2}' Manta20Mb.sort.merge.bed > Manta20Mb.sort.merge.size.bed
bedtools intersect -f 0.8 -wao -a $SigBED_size -b Manta20Mb.sort.merge.size.bed | less -S
bedtools intersect -f 0.8 -wao -a $SigBED_size -b Manta20Mb.sort.merge.size.bed | sort -k1,1n > CNVnatorSig_Manta20Mb.bed
# add CNVnator CNVRs:/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/All.tmp.sort.bed
bedtools merge -c 1 -o count -i /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/All.tmp.sort.bed | sort -k1,1n | awk -v OFS='\t' '{print $0, $3-$2}' > CNVnatorCNVRs.count.size.bed
bedtools intersect -f 0.8 -wao -a CNVnatorSig_Manta20Mb.bed -b CNVnatorCNVRs.count.size.bed > CNVnatorSig_Manta20Mb_CNVnator.bed
### FINAL DATA
cp CNVnatorSig_Manta20Mb_CNVnator.bed /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/Both/
#bedtools intersect -f 0.8 -wao -a $SigBED -b All.filter.tmp.bed | less -S

less -S /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/Population_collapse.bed
#awk -v OFS='\t' '($2 == 54558100)' chr19_54399999_54599999.bed | awk '{ {count = 0}; for (i = 4; i <= NF; i++){ if ($i == 2) {Twocount++}; if (0.5 < $i && $i < 2) {Onecount++}; if ($i < 0.5) {Zerocount++}}; print "Two copies =" Twocount, "One copy =" Onecount, "Zero copy =" Zerocount }'
module load R/3.4.3
module load R_packages/3.4.3
# in R
#cd /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/
#CNVtable <- read.table("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/Population_collapse.bed", header = FALSE, sep = "\t")
CNVtable <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNVbed.rds")
AgeSex <- read.csv(file="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv", header=TRUE, sep=",")
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
rownames(CNVtable) <- paste0(CNVtable$CHROM, ":", CNVtable$Start, "-", CNVtable$End)
# 取biomar与CNV的sample交集: use "%in%":
select.CNVtable <- CNVtable[ , names(CNVtable) %in% rownames(pea_3) ]
# Check size: length(new_output); dim(new_output)
# 取new_testCNV与AgeSex的sample交集: use "%in%"
rownames(AgeSex) <- AgeSex$id #OR rownames(new_AgeSex) <- new_AgeSex[,1]
select.AgeSex <- AgeSex[rownames(AgeSex) %in% colnames(select.CNVtable), ]
# Biomarker has 892 samples, select to 872 samples:
#select.testPea_3 <- testPea_3[ rownames(testPea_3) %in% colnames(select.testCNV), ]
# set the low quality variants to "NA"
select.CNVtable[select.CNVtable<0] <- NA # set negative values to 'NA'
select.CNVtable <- select.CNVtable[rowSums(is.na(select.CNVtable))/length(select.CNVtable)<0.1, ]
dim(select.CNVtable)
# 6 June
library(stringr)
#Coor_CNV <- data.frame( Coor=(rownames(select.CNVtable)) )
library(tidyr)
Coor_CNV <- str_split_fixed(rownames(select.CNVtable), "[:]", 2)
Coor_StartEnd <- str_split_fixed(Coor_CNV[,2], "[-]", 2)
withCoor <- data.frame( Chr=Coor_CNV[,1], Star=Coor_StartEnd[,1], End=Coor_StartEnd[,2], select.CNVtable )
# Save the copy number genotype in the population, did not save the header line: did not save the sample.
write.table(withCoor, "SelectCNV.txt", sep="\t", col.names=FALSE, row.names=FALSE)

# Counting the copy number genotype
SigBED_size=/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed
cp SelectCNV.txt AllCNVnator_CN.bed
# remove "
sed -i "s/\"//g" AllCNVnator_CN.bed
AllCN=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/AllCNVnator_CN.bed
bedtools intersect -wao -a $SigBED_size -b $AllCN | less -S
bedtools intersect -wao -a $SigBED_size -b $AllCN > sig_AllCN.bed
awk -v OFS='\t' '{ {Morecount = 0; Twocount = 0; Onecount = 0; Zerocount = 0}; for (i = 8; i <= NF; i++){ if ($i == 2) {Twocount++}; if (0.5 < $i && $i < 2) {Onecount++}; if ($i < 0.5) {Zerocount++}; if ($i > 2) {Morecount++}}; print $1,$2,$3,$4,$5,$6,$7,$7-$6,"OverTWO:"Morecount, "TWO:"Twocount, "ONE:"Onecount, "Zero:"Zerocount }' sig_AllCN.bed | less -#!/bin/sh
# Function to estimate the copy number genotype in the population
awk -v OFS='\t' '{ {Morecount = 0; Twocount = 0; Onecount = 0; Zerocount = 0}; for (i = 8; i <= NF; i++){ if ($i == 2) {Twocount++}; if (0.5 < $i && $i < 2) {Onecount++}; if ($i < 0.5) {Zerocount++}; if ($i > 2) {Morecount++}}; print $1,$2,$3,$4,$5,$6,$7,$7-$6,Morecount":"Twocount":"Onecount":"Zerocount }' sig_AllCN.bed > sig_AllCN_count.bed
sort -k1,1n  sig_AllCN_count.bed > sig_AllCN_count_sort.bed

#Coor_CNV[,1] <- substring(Coor_CNV[,1], 2)
#write.table(select.CNVtable, "SelectCNV.txt", sep="\t")
#243987    872
library(ggplot2)
All.list <- data.frame()
d1 <- data.frame(CN=unlist(select.CNVtable, use.names = FALSE))
library(ggplot2)
# Basic histogram
ggplot(d1, aes(x=CN)) + geom_histogram()

aspect_ratio <- 1.5
height <- 9
ggsave( "CN.png", height = height, width = height * aspect_ratio)
# set the CN over 10 copies to 10  "
d1.filter <- apply(d1, 1:2, function(x) ifelse(x > 10, 10, x))

saveRDS(d1.filter,"CN.rds")
d1.filter <- readRDS("CN.rds") #  "
d1.filter <- as.data.frame(d1.filter)

#colnames(d1.filter) <-
library(ggplot2)
ggplot(d1.filter, aes(x=CN)) + geom_histogram()
# ggplot(d1.filter, aes(x=CN)) + geom_histogram( binwidth = 0.1)
aspect_ratio <- 1.5
height <- 9
ggsave("CN_10fix.png", height = height , width = height * aspect_ratio) # "

ggplot(d1.filter, aes(x=CN)) + geom_histogram( binwidth = 0.1) +
scale_x_continuous(name = "Copy Number",
                           breaks = seq(0, 10, 1),
                           limits=c(0, 10.5)) +
scale_y_continuous(limits=c(0, 50000000))
ggsave("CN_10fix_5e7.png", height = height , width = height * aspect_ratio)  # "

ggplot(d1.filter, aes(x=CN)) + geom_histogram( binwidth = 0.1) +
scale_x_continuous(name = "Copy Number",
                           breaks = seq(0, 10, 1),
                           limits=c(0, 10.5))
ggsave("CN_10fix_xlab.png", height = height , width = height * aspect_ratio) # "

ggplot(d1.filter, aes(x=CN)) + geom_histogram( binwidth = 0.1) +
scale_x_continuous(name = "Copy Number",
                           breaks = seq(0, 10, 1),
                           limits=c(0, 10.5)) + ylim(0, 2500000)
ggsave("CN_10fixXY_xlab.png", height = height , width = height * aspect_ratio)
