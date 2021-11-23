# Script for summer project:
# Check in the middle of a deletion variants to estimate the number of sample having 2/1/0 deletion for this variant. (number of 0/0 0/1 1/1 genotype indeviduals for the example figures on Chr9)
awk -v OFS='\t' '($2 == 54558100)' chr19_54399999_54599999.bed | awk '{ {count = 0}; for (i = 4; i <= NF; i++){ if ($i == 2) {Twocount++}; if (0.5 < $i && $i < 2) {Onecount++}; if ($i < 0.5) {Zerocount++}}; print "Two copies =" Twocount, "One copy =" Onecount, "Zero copy =" Zerocount }'
#Two copies =713 One copy =291 Zero copy =27

#The CNP-biomarker GLM association p-values for all CNPloci~biomarkers(243,987*438):
/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/T_All_P_only.rds
#All stats results: Estimate/Std.Error/t.value/Pro(>|t|):
/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/WG_P_tables/T_AllPvalues.rds

#Signicifant results
#All results passed the adjusted p-value (0.05/438/245987) threshold are listed in:
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.bed
#17 accociated biomarkers whole genome P-values that extracted for making Manhattan plot using qqman:manhattan function:
less -S  /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/Manha17.csv

#The significant CNVRs and the copy number genotype frequency ($9=Over2:2:1:0):
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/sig_AllCN_count_sort.bed
#Significant CNVRs with copy number genotype:
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/sig_AllCN.bed

#3July:
# add header line for /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/AllCNVnator_CN.bed
# tried bedtools intersect in R?
cd /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/
module load R/3.4.3
module load R_packages/3.4.3
# in R: CNP_filtr.txt is just the CNP table with basic 10% high quality filter
withCoor_fixNames <- data.frame( Chr=Coor_CNV[,1], Star=Coor_StartEnd[,1], End=Coor_StartEnd[,2], select.CNVtable, check.names = FALSE )
write.table(withCoor_fixNames, "CNP_filtr.txt", sep="\t", col.names=TRUE, row.names=FALSE)
#CNP_filtr.txt is SelectCNV.txt in 5_SigSNVRs.sh with header
# in Bash: remove quote(""): sed -i "s/\"//g" /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt
# used vim to add # for the header line
# the above steps make reading the txt file in R return error: more columns than column names
#### Population CNP matrix with header: ####
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt

# TODO: format list of significant CNVRs with associated biomarkers before bedtools intersect:
# raw table: less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.csv
sed 's/,/\t/g' List_significantBiomarkers_withStat.csv > List_significantBiomarkers_withStat.bed
# formatted one, need to reformat(remove some columns): less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.bed

#9July:
cd /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute
#remove the quote("")
sed -i "s/\"//g" /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/List_significantBiomarkers_withStat.bed
#keep 2 3 4 6 8 columns, format as correct bed file
# $4:Estimate/Std.Error/t-value/Pr(>|t|)
cut -d$'\t' -f 2-4,6,8 ../GWAS/List_significantBiomarkers_withStat.bed > ../GWAS/List_significantBiomarkers_withStat_cut.bed
# Add size of the
awk -v OFS='\t' '{print $1, $2, $3, $3-$2, $4, $5}' ../GWAS/List_significantBiomarkers_withStat_cut.bed > ../GWAS/List_significantBiomarkers_withStat_cutSize.bed
module load bioinfo-tools BEDTools/2.27.1
# move the first 4 columns(coor of CNVRs) after $5-7
cut -f 5-7,1-4,8- #! can not change the order
# reorder, for bedtools merging the window result to window copy number genotype:
## for window bigger than 200 bp, they were merged from adjecent windows having same copy number genotype acros the population:
# move the window coordiante to the first 3 columns, before were the significant region coordinate
awk -v OFS='\t' '{print $5,$6,$7,$0}' sig_AllCN.bed > sig_AllCN_reOrder.bed
mkdir ../GWAS/Sig_CNVR_info
# no sample name info, maybe don't need the sample names? header line in:  /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt
# TODO: change ../GWAS/List_significantBiomarkers_withStat_cut.bed to ../GWAS/List_significantBiomarkers_withStat_cutSize.bed
#bedtools intersect -wao -b sig_AllCN_reOrder.bed -a ../GWAS/List_significantBiomarkers_withStat_cut.bed > ../GWAS/Sig_CNVR_info/Intrsect.bed
#bedtools intersect -wao -b CNP_filtr.txt -a ../GWAS/List_significantBiomarkers_withStat_cut.bed | less -S
cd /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute
bedtools intersect -wao -b sig_AllCN_reOrder.bed -a ../GWAS/List_significantBiomarkers_withStat_cutSize.bed > ../GWAS/Sig_CNVR_info/Intrsect.bed
# Intersect result of sig_biomarkers, window and collpased region, Population Copy Number genotype:
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/Sig_CNVR_info/Intrsect.bed

# make more compact, remove some columns
# add header line, in R? by column name function
# Todo: visualiztion
# Create a loop to run through each of the 30 significant CNVRS, creat 30 folder and make histogram and linear regression to the associated biomarkers

# 22 JUL, read Intrsect.bed
# the last 2 columns is the size of the window. they are from two times of bedtools intersect results
