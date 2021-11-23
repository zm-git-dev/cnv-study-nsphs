# reformat the script to debug
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder

DataDir <- file.path("/proj/sens2016007/nobackup/Zhiwei")
CNVtable <- readRDS(file.path(DataDir, "BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNVbed.rds"))
#R > dim(CNVtable)
# [1] 263160   1024
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
rownames(CNVtable) <- paste0(CNVtable$CHROM, ":", CNVtable$Start, "-", CNVtable$End)
select.CNVtable <- CNVtable[ , names(CNVtable) %in% rownames(pea_3) ]

select.CNVtable[select.CNVtable<0] <- NA # set negative values to 'NA'
select.CNVtable <- select.CNVtable[rowSums(is.na(select.CNVtable))/length(select.CNVtable)<0.1, ]
dim(select.CNVtable)
# [1] 243987    872

##new 26 Aug: count the number of regions in the table: Select_CNVtable_withName.bed is generated in next session.
cp /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed ./
# number of line
AllCN="Select_CNVtable_withName.bed"
module load bioinfo-tools BEDTools/2.27.1
awk '{OFS="\t"; print $1, $2, $3}' Select_CNVtable_withName.bed | bedtools merge | wc -l
#29400 number of individual: 872. Population size different for each biomarker association.
## above new 26 Aug

# order by sample name:
select.CNVtable <- select.CNVtable[ ,order(names(select.CNVtable))]

library(stringr)
# Coor_CNV <- data.frame( Coor=(rownames(select.CNVtable)) )

library(tidyr)
Coor_CNV <- str_split_fixed(rownames(select.CNVtable), "[:]", 2)
Coor_StartEnd <- str_split_fixed(Coor_CNV[,2], "[-]", 2)
withCoor <- data.frame( Chr=Coor_CNV[,1], Star=Coor_StartEnd[,1], End=Coor_StartEnd[,2], select.CNVtable, check.names = FALSE )
# Save the copy number genotype in the population, with the sample names:
# selective: select by 10% high quality genotype
write.table(withCoor, "Select_CNVtable_withName.txt", sep="\t", col.names=TRUE, row.names=FALSE)
# in Bash: # format CNV selected matrix:
cp Select_CNVtable_withName.txt Select_CNVtable_withName.bed
# remove "
sed -i "s/\"//g" Select_CNVtable_withName.bed

# CN matrix, whole gnome:
AllCN=/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed
wharfFolder="/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007"
# the 30 sifnificant CNVRs
SigCNRV_size=$wharfFolder/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed
#module load bioinfo-tools BEDTools/2.27.1
#bedtools intersect -wao -a $SigBED_size -b $AllCN | less -S
#bedtools intersect -wao -a $SigBED_size -b $AllCN > sig_AllCN.bed


## in Bash:
cp /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.csv ./
sed 's/,/\t/g' List_significantBiomarkers_withStat.csv > List_significantBiomarkers_withStat.bed
#remove the quote("")
sed -i "s/\"//g" List_significantBiomarkers_withStat.bed
#remove columns for coordinate short name
cut -d$'\t' -f 2-6,8 List_significantBiomarkers_withStat.bed > List_significantBiomarkers_withStat_cut.bed
#add window size:
awk -v OFS='\t' '{print $1, $2, $3, $3-$2, $4, $5, $6}' List_significantBiomarkers_withStat_cut.bed > List_significantBiomarkers_withStat_cutSize.bed
#remove the first line: Chr     Start   End     0       P_value Stat    Biomarker
tail -n +2 List_significantBiomarkers_withStat_cutSize.bed > List_significantBiomarkers_withStat_cutSize.tmp.bed && mv List_significantBiomarkers_withStat_cutSize.tmp.bed List_significantBiomarkers_withStat_cutSize.bed


########BEDtools process: !
# CN matrix, whole gnome:
AllCN="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed"
wharfFolder="/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007"
# the 30 sifnificant CNVRs
SigCNRV_size=$wharfFolder/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed
# the 382 indepedent assiciation:
SigAssociation="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/List_significantBiomarkers_withStat_cutSize.bed"
module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -wao -a $SigAssociation -b $AllCN | wc -l
#382
wc -l $SigAssociation
#382
bedtools intersect -wao -a $SigAssociation -b $AllCN > SigAssociation_CN.bed
# todo:remove the last column
bedtools intersect -wao -a $SigCNRV_size -b SigAssociation_CN.bed > SigCNVR_SigAssociation_CN.bed
# todo: remove the last two columns!

## 26 Aug: lower p-value to include more significant regions:

# continue in Start_changeP_reFormatCNVmatrix.sh
## 26 Aug: extend the window:
wharfFolder="/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007"
# the 30 sifnificant CNVRs:
SigCNRV_size=$wharfFolder/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed
# extend windows by 2000 bp on each side:
awk '{OFS="\t"; print $1, $2-2000, $3+2000}' $SigCNRV_size
