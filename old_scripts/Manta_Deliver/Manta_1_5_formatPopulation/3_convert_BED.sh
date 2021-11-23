#!/bin/bash
module load bioinfo-tools BEDTools/2.27.1
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllBED
mkdir -p $OutPATH
rm -f $OutPATH/*
# save the VCF files directory to one txt file.
#ls -d  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_*/IntermediaOutput/*_Folder/results/variants/diploidSV.vcf.gz > $(dirname $OutPATH)/sample_files_VCF.txt
# use the file in the VCF folder made in 1_moveVCF?
# Since I need to use BEDtools multiple times, using Python will be more complicated?
FILES="/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllVCF/*"
for f in $FILES
do
  BaseName=$(basename ${f%%.*})
  echo Converting "$BaseName" to BED and set Copy Number, seperate regions with overlapping variants.
  zcat $f | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe | \
  awk '(($11=="DUP"||$11=="DEL")&&$12=="PASS"){OFS="\t";
  if($1==$4&&$6>$2&&($6-$2)<=2000000){print $1,$2,$6,$6-$2,$11,$12,substr($15, 1, 3)};
  if($1==$4&&$6<$2&&($3-$5)<=2000000){print $1,$5,$3,$3-$5,$11,$12,substr($15, 1, 3)};
  }' > ${OutPATH}/${BaseName}.PASS.2M.bed
  #remove MT chromosome record so that BEDTools can finish with error: but still other unassembled regions in the bed.
  sed  -i -e '/MT/d' ${OutPATH}/${BaseName}.PASS.2M.bed
  #Unable to tell the sex just using the X and Y chromosome record, #TODO: special format for sex chromosome
  sort -k1,1 -k2,2n ${OutPATH}/${BaseName}.PASS.2M.bed | \
   bedtools merge -i stdin -c 1 -o count | awk '$4 == 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | \
   bedtools intersect -wao -a stdin -b ${OutPATH}/${BaseName}.PASS.2M.bed | \
   awk '{OFS="\t";
   if($10=="DEL"&&($12 ~ /^1\/1/) ){print $6,$7,$8,$9,$10,$11,$12,0};
   if($10=="DEL"&&($12 ~ /^0\/1/) ){print $6,$7,$8,$9,$10,$11,$12,1};
   if($10=="DUP"&&($12 ~ /^1\/1/) ){print $6,$7,$8,$9,$10,$11,$12,4};
   if($10=="DUP"&&($12 ~ /^0\/1/) ){print $6,$7,$8,$9,$10,$11,$12,3};
   }' > ${OutPATH}/${BaseName}.PASS.2M.sort.rmOverlap.bed
   # copy number variable regions with over one variants record,
   #TODO: need to improve: regions with two DEL and DUP set as 'NA'
  sort -k1,1 -k2,2n ${OutPATH}/${BaseName}.PASS.2M.bed | \
   bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,"#Overlap:"$4,".",".","NA"}' > ${OutPATH}/${BaseName}.PASS.2M.sort.OVERLAP.bed
  cat ${OutPATH}/${BaseName}.PASS.2M.sort.rmOverlap.bed ${OutPATH}/${BaseName}.PASS.2M.sort.OVERLAP.bed > ${OutPATH}/${BaseName}.PASS.2M.sort.General.bed
done
# for the downstream analysis combine reOverlap and OVERLAP data?

# sort -k1,1 -k2,2n /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllBED/KA09-155.PASS.2M.bed | \
# bedtools merge -i stdin -c 1 -o count | awk '$4 > 1 {OFS="\t"; print $1,$2,$3,$3-$2,$4}' | \
# bedtools intersect -wao -a stdin -b /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllBED/KA09-155.PASS.2M.bed | less -S
#sbatch -A sens2016007 -p core -n 4 -t 2:00:00 -J
