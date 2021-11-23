#load surviror and vcftools
module load bioinfo-tools vcftools/0.1.15
module load bioinfo-tools SURVIVOR/1.0.3

#awk -F '\t' '{if($0 ~ /\#\) print; else if($3 ~ /MantaDEL/) print}'
# $10 start wit 1, homozygous DEL DUP.
awk -F '\t' '{if($0 ~ /\#/) print; if((($3 ~ /^MantaDEL/)||($3 ~ /^MantaDUP/))&&($10 ~ /^1\/1/)) print}' P2656_301.vcf  | less -S

# use one of the PacBio seq individual to get the numbers of overlapping regions:
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/TwoPacBioSamples/KA06_0017_P2656_342_Folder
cd test_overlap_forPopulation
module load bioinfo-tools BEDTools/2.27.1
# regions plan to exclude:

sort -k1,1 -k2,2n ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge | bedtools intersect -wao  -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk -F '\t' '{print $1,$2,$3}' | uniq -c | awk '$1 > 1 {OFS="\t"; print $2,$3,$4,$1}'
# okey to remove all the multiple observations regions?
sort -k1,1 -k2,2n ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge | bedtools intersect -wao  -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk -F '\t' '{print $1,$2,$3}' | uniq -c | awk '$1 > 1 {OFS="\t"; print $2,$3,$4,$4-$3,$1}' | bedtools intersect -wao -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | less -S
# change to numeric chr name for the repeat.bed
#sed 's/^chr\|%$//g' file
#/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/RepeatMaster/hg19_repeat.bed
wc -l ../Manta_P2656_342.PASS_Ploidy_Size2M.bed
# 4656
# exclude the regions with multiple variants
sort -k1,1 -k2,2n ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge | bedtools intersect -wao  -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk -F '\t' '{print $1,$2,$3}' | uniq -c | awk '$1 == 1 {OFS="\t"; print $2,$3,$4,$4-$3}' | bedtools intersect -wao -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk '{print $10}' | sort | uniq -c
#   4393 PASS
# use bedtools merge -n to get the number of overlap instead:
sort -k1,1 -k2,2n ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge -i stdin -c 1 -o count | awk '$4 == 1 {OFS="\t"; print $1,$2,$3,$3-$2}' | bedtools intersect -wao -a stdin -b ../Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk '{print $10}' | sort | uniq -c
#   4393 PASS
#NO! try a Python function to format the output bed file in chr, start, end, CNestimation
# save as above, convert the variantType+GT to copy number extimation, add number of variants in the region for testing:
sort -k1,1 -k2,2n $eachBED | head -n10 | \
 bedtools merge -i stdin -c 1 -o count | awk '$4 == 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | \
 bedtools intersect -wao -a stdin -b $eachBED

eachBED="../Manta_P2656_342.PASS_Ploidy_Size2M.bed"
sort -k1,1 -k2,2n $eachBED | head -n10 | \
 bedtools merge -i stdin -c 1 -o count | awk '$4 == 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | \
 bedtools intersect -wao -a stdin -b $eachBED | \
 awk '{OFS="\t";
 if($10=="DEL"&&($13 ~ /^1\/1/) ){print $6,$7,$8,$9,$10,substr($13, 1, 3),0,$11};
 if($10=="DEL"&&($13 ~ /^0\/1/) ){print $6,$7,$8,$9,$10,substr($13, 1, 3),1,$11};
 if($10=="DUP"&&($13 ~ /^1\/1/) ){print $6,$7,$8,$9,$10,substr($13, 1, 3),4,$11};
 if($10=="DUP"&&($13 ~ /^0\/1/) ){print $6,$7,$8,$9,$10,substr($13, 1, 3),3,$11};
}'

# multipe variants regions marked as NA:, also region with two dups
OutPATH=/proj/sens2016007/nobackup/Zhiwei/SummerProject/TwoPacBioSamples/KA06_0017_P2656_342_Folder
sort -k1,1 -k2,2n ${OutPATH}/Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | bedtools intersect -wao -a stdin -b ${OutPATH}/Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk '{print $11}' | sort | uniq -c
sort -k1,1 -k2,2n ${OutPATH}/Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | bedtools intersect -wao -a stdin -b ${OutPATH}/Manta_P2656_342.PASS_Ploidy_Size2M.bed > ${OutPATH}/Manta_P2656_342.PASS_Ploidy_Size2M_OverlapDetail.bed
# solve overlapping region:
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/TwoPacBioSamples/KA06_0017_P2656_342_Folder
sort  -k1,1 -k2,2n Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' \
| bedtools intersect -wao -a stdin -b Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$10}' | uniq -c | awk '($1 == $6) {OFS="\t";print $1,$2}' |
# not working! All, overlap with general region:
sort -k1,1 -k2,2n Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | bedtools intersect -wao -a stdin -b Manta_P2656_342.PASS_Ploidy_Size2M.bed | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$10}' | uniq -c | awk '($1 == $6) {OFS="\t"; print $2,$3,$4}' | bedtools subtract -a stdin -b Manta_P2656_342.PASS_Ploidy_Size2M.bed
#TODO: get the overlaping and non-overlapping region, mark as homozygous and hetezygous.
#211 PASS
# 52 Ploidy
# the converted bed file columns: Region_chr, Region_start, Region_end, Event:NA, GT:
