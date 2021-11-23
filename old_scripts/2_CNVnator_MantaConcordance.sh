###Filter, summarize and compare population-wide CNVnator, Manta and dbSV(A Copy Number Variation Map of the Human Genome (Nature Reviews Genetics, 2015) Gain+Loss hg19) results.
# Compare intersect, in BASH
OutPATH=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest # MORE Stringent filter quiteria $9:q0<0.5 and $5:e-val1(t-test)<0.05/1047
mkdir -p $OutPATH
# list of the raw data folder:
VariantTXT=$(ls -d /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/3_Variant_output/*)
for CNVnatorCall in ${VariantTXT}
do
  BaseName=${CNVnatorCall##*/}
  awk '{if ($9>=0&&$9<=0.5&&$5<=0.00005) print $0}' $CNVnatorCall | awk -v OFS="\t" '{split($2, a, /:|-/); print a[1], a[2]-1, a[3], $4*2, $5, $6, $9, $3}' \
  > ${OutPATH}/${BaseName%.txt}_selected.bed
done
# Count variant after qo+t_test filter, mean value of observation for each sample
awk '($4>=2){DUPsum+=1}($4<2){DELsum+=1} END {print "DEL mean =", DELsum/1047, "DUP mean =", DUPsum/1047, "Average =" (DELsum+DUPsum)/1047}' /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/*
#DEL mean = 744.171 DUP mean = 640.695 Average =1384.87
# CNVnator stat to be compared with Manta:
OutPATH=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder
mkdir -p $OutPATH
cd $OutPATH
module load bioinfo-tools BEDTools/2.27.1
# get all variant from bed files in the OptBinSize_q0_Ttest folder
cat ../*_selected.bed > All.tmp.bed
sort -k1,1 -k2,2n All.tmp.bed > All.tmp.sort.bed
# Split by DUP DEL for comparing the CNVRs seperately for DUP and DEL.
awk -v OFS="\t" '($4 > 2){print $1, $2, $3, $4, $8}' All.tmp.sort.bed > All.tmp.sort.DUP.bed
awk  -v OFS="\t" '($4 <= 2){print $1, $2, $3, $4, $8}' All.tmp.sort.bed > All.tmp.sort.DEL.bed
awk -v OFS="\t" '{print $1, $2, $3, $4, $8}' All.tmp.sort.bed > All.tmp.sort.short.bed
#TODO: plot them using R DUP+DEL for CNVnator only, All.tmp.sort.bed for comparing CNVnator and Manta.
## NEW 21April, merge to reduce repeat, add size of the variant in $5, can be used by the R scripts used before for All.tmp.sort.DUP/DEL.bed
$ module load bioinfo-tools BEDTools/2.27.1
$ bedtools merge -i All.tmp.sort.short.bed -c 1 -o count | awk -v OFS="\t" '{print $0, $3-$2}' > All.tmp.sort.CNP.mergeCount.bed
$ bedtools merge -i All.tmp.sort.DUP.bed -c 1 -o count | awk -v OFS="\t" '{print $0, $3-$2}' > All.tmp.sort.DUP.mergeCount.bed
$ bedtools merge -i All.tmp.sort.DEL.bed -c 1 -o count | awk -v OFS="\t" '{print $0, $3-$2}' > All.tmp.sort.DEL.mergeCount.bed
# compare with know SVdb in /proj/sens2016007/nobackup/Zhiwei/SVdb
Gain_db=/proj/sens2016007/nobackup/Zhiwei/SVdb/Inclusive.Gain.bed
Loss_db=/proj/sens2016007/nobackup/Zhiwei/SVdb/Inclusive.Loss.bed
bedtools intersect -v -a All.tmp.sort.DUP.mergeCount.bed -b $Gain_db
# -v reverse, count the variants that are not in -b file
wc -l $Gain_db
#1983 #1983-1829=154; 1983-1603=380
#bedtools intersect -v -r -f 0.5  -b All.tmp.sort.DUP.mergeCount.bed -a $Gain_db | wc -l #1829
#bedtools intersect -v -b All.tmp.sort.DUP.mergeCount.bed -a $Gain_db | wc -l #1603
wc -l All.tmp.sort.DUP.mergeCount.bed #2805   2651/2805=0,9450980392       # 2805-2651=154; 2805-2153=652==> different from above, should use minimum 50% overlap?
bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DUP.mergeCount.bed -b $Gain_db | wc -l #2651
wc -l All.tmp.sort.DEL.mergeCount.bed #8579 6578/8579=0,7667560322
bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DEL.mergeCount.bed -b $Loss_db | wc -l #6578
#bedtools intersect -v -a All.tmp.sort.DUP.mergeCount.bed -b $Gain_db | wc -l #2153

# Compare all variants
wc -l All.tmp.sort.DEL.bed #779147 531380/779147 = 0,6820022409
bedtools intersect -v -r -f 0.5 -a All.tmp.sort.DEL.bed -b $Loss_db | wc -l #531380
wc -l All.tmp.sort.DUP.bed #670808 649059/670808 = 0,967577906
bedtools intersect -v -r -f 0.5 -a All.tmp.sort.DUP.bed -b $Gain_db | wc -l #649059


##Compare CNVnator merge and Manta merge: check 5_SigCNVRs.sh for Manta1.5 formatting

# Manta filter by PASS only, Update 1 May: TODO: just the result of Manta1.5, it should be finished, check 5_SigCNVRs.sh for scripts used to format Manta1.5 results.
# OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/MantaMain/results_allBED_onlyPASS
# mkdir -p $OutPATH
# VariantGZ=$(ls -d  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/MantaMain/Manta_xa*/IntermediaOutput/*_Folder/results/variants/diploidSV.vcf.gz)
# for Call in ${VariantGZ}
# do
#   BaseName=$(basename $(echo ${Call} | sed 's/.*\(IntermediaOutput\)/\1/g' | cut -d'/' -f-2) )
#   echo Saving to ${OutPATH}/${BaseName%_Folder}.PASS.bed
#   zcat $Call | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe | awk '(($11=="DUP"||$11=="DEL")&&$12=="PASS"&&$1==$4&&$6>$2){OFS="\t";print $1,$2,$6,$6-$2,$0}(($11=="DUP"||$11=="DEL")&&$12=="PASS"&&$1==$4&&$6<$2){OFS="\t";print $1,$5,$3,$3-$5,$0}' > ${OutPATH}/${BaseName%_Folder}.PASS.bed
# done
# # Count variant after PASS filter
# awk '($15=="DUP"){DUPsum+=1}($15=="DEL"){DELsum+=1} END {print "DEL mean =", DELsum/1047, "DUP mean =", DUPsum/1047, "Average =" (DELsum+DUPsum)/1047}' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/MantaMain/results_allBED_onlyPASS/*
# #DEL mean = 4600.35 DUP mean = 481.007 Average =5081.36
# OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/MantaMain/results_allBED_onlyPASS/StatFolder
# mkdir -p $OutPATH
# cd $OutPATH
# module load bioinfo-tools BEDTools/2.27.1
# #cat ../*PASS.bed > All.tmp.bed
# awk -v OFS="\t" '{print $1, $2, $3, $15, $4}' ../*PASS.bed > All.tmp.bed
# sort -k1,1 -k2,2n All.tmp.bed > All.tmp.sort.bed
# awk -v OFS="\t" '($4 == "DUP"){print $0}' All.tmp.sort.bed > All.tmp.sort.DUP.bed
# awk -v OFS="\t" '($4 == "DEL"){print $0}' All.tmp.sort.bed > All.tmp.sort.DEL.bed
# # Chech the variatns size distribution, there are some BIG variants (1000Mb) reported by MANTA!
# awk '{print $5}' All.tmp.sort.bed | sort -k1,1n | uniq -c | less -S
# #TODO: plot them using R DUP+DEL for Manta only, All.tmp.sort.bed for comparing CNVnator and Manta.
# #21April new, merge:
# $ module load bioinfo-tools BEDTools/2.27.1
# # remove MT seq:
# sed  -i -e '/MT/d' All.tmp.sort.DUP.bed
# bedtools merge -i All.tmp.sort.DUP.bed -c 1 -o count | awk -v OFS="\t" '{print $0, $3-$2}' > All.tmp.sort.DUP.mergeCount.bed
# bedtools merge -i All.tmp.sort.DEL.bed -c 1 -o count | awk -v OFS="\t" '{print $0, $3-$2}' > All.tmp.sort.DEL.mergeCount.bed
# natorDUP=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder/All.tmp.sort.DUP.mergeCount.bed
# natorDEL=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder/All.tmp.sort.DEL.mergeCount.bed
# wc -l $natorDUP #2805-2736=69 2736/2805=0,9754010695
# wc -l All.tmp.sort.DUP.mergeCount.bed #1310-1241=69 1241/1310=0,9473282443
# bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DUP.mergeCount.bed -b $natorDUP | wc -l #1241
# bedtools intersect -v -r -f 0.5  -b All.tmp.sort.DUP.mergeCount.bed -a $natorDUP | wc -l #2736
# wc -l $natorDEL #8579-8140=439 8140/8579 =0,9488285348
# wc -l All.tmp.sort.DEL.mergeCount.bed # 7564-7125=439 7125/7564=0,9419619249
# bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DEL.mergeCount.bed -b $natorDEL | wc -l # 7125
# bedtools intersect -v -r -f 0.5  -b All.tmp.sort.DEL.mergeCount.bed -a $natorDEL | wc -l # 8140
#
# # Reduce the reciprocal threshold to 20%
# wc -l $natorDUP #2805-2714= 2714/2805=0,9675579323
# wc -l All.tmp.sort.DUP.mergeCount.bed #1310-1225=85 1225/1310=0,9351145038
# bedtools intersect -v -r -f 0.2  -a All.tmp.sort.DUP.mergeCount.bed -b $natorDUP | wc -l #1225
# bedtools intersect -v -r -f 0.2  -b All.tmp.sort.DUP.mergeCount.bed -a $natorDUP | wc -l #2714
# wc -l $natorDEL #8579-8013= 8013/8579 =0,9340249446
# wc -l All.tmp.sort.DEL.mergeCount.bed # 7564-6995= 6995/7564=0,9247752512
# bedtools intersect -v -r -f 0.2  -a All.tmp.sort.DEL.mergeCount.bed -b $natorDEL | wc -l # 6995
# bedtools intersect -v -r -f 0.2  -b All.tmp.sort.DEL.mergeCount.bed -a $natorDEL | wc -l # 8013
#
#
# # Compare to database:
# Gain_db=/proj/sens2016007/nobackup/Zhiwei/SVdb/Inclusive.Gain.bed
# Loss_db=/proj/sens2016007/nobackup/Zhiwei/SVdb/Inclusive.Loss.bed
# wc -l $Gain_db # 1983
# wc -l All.tmp.sort.DUP.mergeCount.bed # 1310  1288/1310=0,9832061069
# bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DUP.mergeCount.bed -b $Gain_db | wc -l #1288
# wc -l All.tmp.sort.DEL.mergeCount.bed # 7564 6647/7564=0,8787678477
# bedtools intersect -v -r -f 0.5  -a All.tmp.sort.DEL.mergeCount.bed -b $Loss_db | wc -l # 6647
# #bedtools intersect -v -a All.tmp.sort.DUP.mergeCount.bed -b $Gain_db | wc -l # 1156
# # interset the all variants All.tmp.sort.DEL.bed All.tmp.sort.DUP.bed
# wc -l All.tmp.sort.DEL.bed #4816566 3717413/4816566 = 0,771797376
# bedtools intersect -v -r -f 0.5 -a All.tmp.sort.DEL.bed -b $Loss_db | wc -l #3717413
# wc -l All.tmp.sort.DUP.bed #503611 481501/503611 = 0,956097067
# bedtools intersect -v -r -f 0.5 -a All.tmp.sort.DUP.bed -b $Gain_db | wc -l #481501
#
# #Both all variants compare:
# mantaDUP=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/MantaTables/All.tmp.sort.DUP.bed #503611 503611-470683=32 928
# mantaDEL=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/MantaTables/All.tmp.sort.DEL.bed #4816566
# natorDUP=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/All.tmp.sort.DUP.bed #670808 670808-644172=26 636
# natorDEL=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/All.tmp.sort.DEL.bed #779147
# bedtools intersect -v -r -f 0.5 -a $mantaDUP -b $natorDUP | wc -l #470683
# bedtools intersect -v -r -f 0.5 -b $mantaDUP -a $natorDUP | wc -l #644172
# bedtools intersect -v -r -f 0.5 -a $mantaDEL -b $natorDEL | wc -l #4190595
# bedtools intersect -v -r -f 0.5 -b $mantaDEL -a $natorDEL | wc -l #373316
MantaAll="/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/StatFinalPre/Manta20Mb.sort.merge.size.bed"
CNVnatorAll="/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder/All.tmp.sort.CNP.mergeCount.bed"


##################
# CNVnator and Manta/1.5 size distribution:
# $ cd /proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/Both/
module load R/3.4.3
module load R_packages/3.4.3
# in R
rm(list = ls())
natorDUPbed <- read.table("/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder/All.tmp.sort.DUP.mergeCount.bed", header = FALSE, sep = "\t")
natorDELbed <- read.table("/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0_Ttest/StatFolder/All.tmp.sort.DEL.mergeCount.bed", header = FALSE, sep = "\t")
mantaDUPbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/All.tmp.sort.DUP.mergeCount.bed", header = FALSE, sep = "\t")
mantaDELbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/results_allBED_onlyPASS/StatFolder/All.tmp.sort.DEL.mergeCount.bed", header = FALSE, sep = "\t")
# Save into one dataframe(rbind) add a column marking the 4 groups:
hisSizeData <- rbind(data.frame( Type=rep("CNVnator(DUP)",nrow(natorDUPbed)), Size=natorDUPbed[,5]), data.frame(Type=rep("CNVnator(DEL)",nrow(natorDELbed)), Size=natorDELbed[,5]), data.frame( Type=rep("Manta(DUP)",nrow(mantaDUPbed)), Size=mantaDUPbed[,5]), data.frame( Type=rep("Manta(DEL)",nrow(mantaDELbed)), Size=mantaDELbed[,5]))
# Filter into two groups, smaller than 1kb or over:
Fit15Kb <- hisSizeData
for(i in 1:nrow(Fit15Kb)) {
  if(Fit15Kb[i,2]>15000) {
    Fit15Kb[i,2] <- 15000
  }
}
ggplot(data= Fit15Kb, aes(x= Size, fill=Type)) +
  geom_histogram() +
  theme(plot.title = element_text(size=22), axis.text = element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12)) +
  ggtitle("Distribution of collapsed CNVRs size") +
  xlab("Size bp") +
  ylab("Count")
aspect_ratio <- 1.5
height <- 9
ggsave("mBothggplot.png", height = height , width = height * aspect_ratio)

ggplot(data= Fit15Kb, aes(x= Size, color=Type)) +
  geom_histogram(fill="white") +
  theme(plot.title = element_text(size=22), axis.text = element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=12)) +
  ggtitle("Distribution of collapsed CNVRs size") +
  xlab("Size bp")
aspect_ratio <- 1.5
height <- 9
ggsave("mBothggplotWhite.png", height = height , width = height * aspect_ratio)

ggplot(data= Fit15Kb, aes(x= Size, color=Type)) +
  geom_density()
ggsave("mBothggplotWithinAllDen.png")

########Testing scripts#######
library("ggplot2")
ggplot(data= hisSizeData[hisSizeData[,2]<=1000,], aes(x= Size, color=Type)) +
  geom_histogram(fill="white")
ggsave("mBothggplotWithin1Kb.png")

Fit16Kb <- hisSizeData
for(i in 1:nrow(Fit16Kb)) {
  if(Fit16Kb[i,2]>16000) {
    Fit16Kb[i,2] <- 16000
  }
}
ggplot(data= Fit16Kb, aes(x= Size, color=Type)) +
  geom_density()
ggsave("mBothggplotWithinAllDen.png")

ggplot(data= Fit16Kb[Fit16Kb[,2]>=1000,], aes(x= Size, color=Type)) +
  geom_histogram(fill="white")
ggsave("mBothggplotWithout1Kb.png")

ggplot(data= hisSizeData, aes(x= Size, color=Type)) +
  geom_histogram(fill="white", binwidth=500)
ggsave("mBothggplotTest.png")
aspect_ratio <- 2
height <- 8
ggsave("mBothggplot.png", height = 8 , width = 8 * aspect_ratio)
