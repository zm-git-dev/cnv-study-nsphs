#module load R/3.4.3
#module load R_packages/3.4.3
cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder
#in R
CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_population_bin10_rowname.rds")
dim(CNVbed)
# [1] 124526   1024
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
CNVtable <- CNVbed
rownames(CNVtable) <- paste0(CNVtable$Chr, ":", CNVtable$Start, "-", CNVtable$End)
select.CNVtable <- CNVtable[ , names(CNVtable) %in% rownames(pea_3) ]
dim(select.CNVtable)
# [1] 124526    872
#select.CNVtable[select.CNVtable<0] <- NA # set negative values to 'NA'
select.CNVtable <- select.CNVtable[rowSums(is.na(select.CNVtable))/length(select.CNVtable)<0.1, ]
dim(select.CNVtable)
# [1] 115736    872
Coor_CNV <- data.frame( Coor=(rownames(select.CNVtable)) )
library(tidyr)
library(stringr)
Coor_CNV <- str_split_fixed(rownames(select.CNVtable), "[:]", 2)
Coor_StartEnd <- str_split_fixed(Coor_CNV[,2], "[-]", 2)
withCoor <- data.frame( Chr=Coor_CNV[,1], Star=Coor_StartEnd[,1], End=Coor_StartEnd[,2], select.CNVtable, check.names = FALSE )
# Save the copy number genotype in the population, with the sample names:
# selective: select by 10% high quality genotype
write.table(withCoor, "Select_CNVtable_withName_Manta.txt", sep="\t", col.names=TRUE, row.names=FALSE)
#in BASH:
module load bioinfo-tools BEDTools/2.27.1
#awk '{OFS="\t"; print $1, $2, $3}' Select_CNVtable_withName_Manta.txt | bedtools merge | wc -l
#remove quptes!
sed s/\"//g Select_CNVtable_withName_Manta.txt > Select_CNVtable_withName_Manta.bed
awk '{OFS="\t"; print $1, $2, $3}' Select_CNVtable_withName_Manta.bed | bedtools merge | wc -l
#30682! used to adjusted multiple testing for Manta
