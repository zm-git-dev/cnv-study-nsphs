cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder


#define #CNPs for multiple test!
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

##new 26 Aug: count the number of regions in the table:
cp /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed ./
# number of line
AllCN="Select_CNVtable_withName.bed"
module load bioinfo-tools BEDTools/2.27.1
awk '{OFS="\t"; print $1, $2, $3}' Select_CNVtable_withName.bed | bedtools merge | wc -l
#29400 number of individual: 872. Population size different for each biomarker association.
# ## above new 26 Aug
#
# # order by sample name:
# select.CNVtable <- select.CNVtable[ ,order(names(select.CNVtable))]
#
# library(stringr)
# # Coor_CNV <- data.frame( Coor=(rownames(select.CNVtable)) )
#
# library(tidyr)
# Coor_CNV <- str_split_fixed(rownames(select.CNVtable), "[:]", 2)
# Coor_StartEnd <- str_split_fixed(Coor_CNV[,2], "[-]", 2)
# withCoor <- data.frame( Chr=Coor_CNV[,1], Star=Coor_StartEnd[,1], End=Coor_StartEnd[,2], select.CNVtable, check.names = FALSE )
# # Save the copy number genotype in the population, with the sample names:
# # selective: select by 10% high quality genotype
# write.table(withCoor, "Select_CNVtable_withName.txt", sep="\t", col.names=TRUE, row.names=FALSE)
# # in Bash: # format CNV selected matrix:
# cp Select_CNVtable_withName.txt Select_CNVtable_withName.bed
# # remove "
# sed -i "s/\"//g" Select_CNVtable_withName.bed

# CN matrix, whole gnome:
AllCN=/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed

# the 39 sifnificant CNVRs
SigCNRV_size="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
#module load bioinfo-tools BEDTools/2.27.1
#bedtools intersect -wao -a $SigBED_size -b $AllCN | less -S
#bedtools intersect -wao -a $SigBED_size -b $AllCN > sig_AllCN.bed


## in Bash:
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
#remove the quote("")
sed "s/\"//g" tmp_format_List_significantBiomarkers_withStat_changePvalue.bed > tmp_format_List_significantBiomarkers_withStat_noQuote_changePvalue.bed
#remove columns for coordinate short name
cut -d$'\t' -f 1-5,7 tmp_format_List_significantBiomarkers_withStat_noQuote_changePvalue.bed > tmp_format_List_significantBiomarkers_withStat_noQuote_cut_changePvalue.bed
rm tmp_format_List_significantBiomarkers_withStat_noQuote_changePvalue.bed
#add window size:
awk -v OFS='\t' '{print $1, $2, $3, $3-$2, $4, $5, $6}' tmp_format_List_significantBiomarkers_withStat_noQuote_cut_changePvalue.bed | sort -k1,1 -k2,2n > List_significantBiomarkers_withStat_cutSize_changePvalue.bed

########BEDtools process: !
# CN matrix, whole gnome:
AllCN="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed"
# the 39 sifnificant CNVRs
SigCNRV_size="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
# the 458 assiciations:
SigAssociation="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/List_significantBiomarkers_withStat_cutSize_changePvalue.bed"
module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -wao -a $SigAssociation -b $AllCN | wc -l
#458
wc -l $SigAssociation
#458
bedtools intersect -wao -a $SigAssociation -b $AllCN > SigAssociation_CN_changeP.bed
# todo:remove the last column
bedtools intersect -wao -a $SigCNRV_size -b SigAssociation_CN_changeP.bed > SigCNVR_SigAssociation_CN_changeP.bed
# todo: remove the last two columns! in R
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/Figures_changeP
####formatting and plot see re_BoxPlot_CNV.r !


##
## 2Sep: lower p-value to include more significant regions:
# TODO: extend the CNVR regions:
## 26 Aug: extend the window:
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder
mkdir ExtendWindow
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow
# CN matrix, whole gnome:
AllCN="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed"
# the 39 sifnificant CNVRs
SigCNRV_size="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
# extend windows by 2000 bp on each side:?, to do with all regions or just the
cat Chr5_70303300.bed
#5       70303300        70312700        9400
awk '{OFS="\t"; print $1, $2-14000, $3+2600, $4+16600}' Chr5_70303300.bed > Chr5_70303300_17kExtend.bed
#module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -wao -a Chr5_70303300_17kExtend.bed -b $AllCN > AllCN_Chr5_70303300_17kExtend.bed
#less -S AllCN_Chr5_70303300_17kExtend.bed
sort -nk6,6 -nk7,7 AllCN_Chr5_70303300_17kExtend.bed > AllCN_Chr5_70303300_17kExtend.sort.bed
# send to R for making plot
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
#in R
list.files()
mybed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow/AllCN_Chr5_70303300_17kExtend.sort.bed", header = FALSE, sep="\t")
# remove the laste column, the size of the overlap
# mybedRm <- subset(mybed, select = -c(880))
# pdf("CNVR_Chr5_70303300.pdf")
# #"
# x <- seq(mybed[1,6],mybed[nrow(mybed),7], length = 300)
# y <- seq(-5, 5, length = 300)
# plot(x,y,type="n",xlim=c(mybed[1,6],mybed[nrow(mybed),7]),ylim=c(0,6), xlab = "Chr5:70290700-70315300,24,600bp", ylab = "Copy number")
# for(j in c(1:nrow(mybed))){
#     for (i in 8:(ncol(mybed)-1)) {
#         # RD <- mybed[i]
#         #str(RD)
#         # segments(x0 = t(startPos), y0 = t(RD), x1 = t(endPos), y1 = t(RD), col=col[i])
#         #print(paste0("row and col factor: ", j, " col:", i))
#         #if  (j > 1) {
#           #segments(x0 = mybed[j-1,7], x1 = mybed[j,6],y0 = mybed[j-1,i], y1 = mybed[j,i] )
#         #}
#         segments(x0 = mybed[j,6], x1 = mybed[j,7],y0 = mybed[j,i], y1 = mybed[j,i] )
#     }
# }
# dev.off()
# in ggplot geom_step:
library(ggplot2)
mybed[is.na(mybed)] <- 2
b <- ggplot(mybed)
for (i in 8:(ncol(mybed)-1)) {
  b <- b + geom_step(data=mybed,mapping=aes_string(x=mybed[,6], y=mybed[,i]))
}
ggsave("CNVR_Chr5_70303300_ggplot.pdf")

b <- ggplot(mybed)
for (i in 8:30) {
  b + geom_step(data=mybed,mapping=aes_string(x=mybed[,6], y=mybed[,i]))
}
ggsave("CNVR_Chr5_70303300_ggplot.pdf")

pdf("CNVR_Chr5_70303300_ggplot.pdf") #"
print(b)
dev.off()
b <- ggplot(mybed) + geom_step(data=mybed, mapping=aes_string(x=mybed[,6], y=mybed[,i]))
b + geom_step(data=mybed, mapping=aes_string(x=mybed[,6], y=mybed[,i+1]))

library(tidyr)
#remove the intersect size info:
mybedRm <- subset(mybed, select = -c(880))
long_mybed <- mybedRm %>% gather(Sample, CN, 8:ncol(mybedRm))

#ggplot(long_mybed)

b <- ggplot(long_mybed,aes(V6,CN))
b + geom_step(aes(group = Sample))
ggsave("CNVR_Chr5_70303300_ggplot.pdf") #"
##GGplot step function
b <- ggplot(long_mybed)
b <- b + geom_step(aes(x=V6,y=CN, group=Sample ))
b + labs(x = "Position", y = "Copy Number", title = "Step Plot")
ggsave("CNVR_Chr5_70303300_ggplot.pdf")


# ggplot(long_mybed, aes(V6,CN, group=Sample, color=Sample)) + geom_step()
#
# #"
# ###new extend the region:
# # in Bash:
# cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow
# module load bioinfo-tools BEDTools/2.27.1
# # CN matrix, whole gnome:
# AllCN="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed"
# awk '{OFS="\t"; print $1, $2-14000, $3+5000, $4+19000}' Chr5_70303300.bed > Chr5_70303300_19kExtend.bed
# #module load bioinfo-tools BEDTools/2.27.1
# bedtools intersect -wao -a Chr5_70303300_19kExtend.bed -b $AllCN > AllCN_Chr5_70303300_19kExtend.bed
# #less -S AllCN_Chr5_70303300_17kExtend.bed
# sort -nk6,6 -nk7,7 AllCN_Chr5_70303300_19kExtend.bed > AllCN_Chr5_70303300_19kExtend.sort.bed
# # in R
# mybed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow/AllCN_Chr5_70303300_19kExtend.sort.bed", header = FALSE, sep="\t")
# mybed[is.na(mybed)] <- 2
# library(tidyr)
# #remove the intersect size info:
# mybedRm <- subset(mybed, select = -c(880))
# long_mybed <- mybedRm %>% gather(Sample, CN, 8:ncol(mybedRm))
# library(ggplot2)
# b <- ggplot(long_mybed)
# b <- b + geom_step(aes(x=V6,y=CN, group=Sample ))
# b + labs(x = "Chr5:70290700-70318300,27,600bp", y = "Copy Number", title = "Copy Number distribution in the Population")
# ggsave("CNVR_Chr5_70303300_ggplot.pdf")
#
# # line plot, not good
# long_mybed_over2 <- long_mybed[long_mybed[,9]!=2,]
# h <- ggplot(long_mybed_over2)
# h <- h + geom_step(aes(x=V6,y=CN, group=Sample ))
# h + labs(x = "Chr5:70290700-70318300,27,600bp", y = "Copy Number", title = "Copy Number distribution in the Population")
# ggsave("CNVR_Chr5_70303300_ggplot_line.pdf")
#
# # also the line plot: #"
# pdf("CNVR_Chr5_70303300.pdf")
# x <- seq(mybed[1,6],mybed[nrow(mybed),7], length = 300)
# y <- seq(-5, 5, length = 300)
# plot(x,y,type="n",xlim=c(mybed[1,6],mybed[nrow(mybed),7]),ylim=c(0,6), xlab = "Chr5:70290700-70318300,27,600bp", ylab = "Copy number", main ="Copy Number distribution in the Population")
# for(j in c(1:nrow(mybed))){
#     for (i in 8:(ncol(mybed)-1)) {
#         # RD <- mybed[i]
#         #str(RD)
#         # segments(x0 = t(startPos), y0 = t(RD), x1 = t(endPos), y1 = t(RD), col=col[i])
#         #print(paste0("row and col factor: ", j, " col:", i))
#         #if  (j > 1) {
#           #segments(x0 = mybed[j-1,7], x1 = mybed[j,6],y0 = mybed[j-1,i], y1 = mybed[j,i] )
#         #}
#         segments(x0 = mybed[j,6], x1 = mybed[j,7],y0 = mybed[j,i], y1 = mybed[j,i] )
#     }
# }
# dev.off()
