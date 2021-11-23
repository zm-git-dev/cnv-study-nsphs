#!/bin/bash -l

#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 12:45:00
#SBATCH -J BinAll_30bp
#SBATCH --mail-type=ALL
#SBATCH --mail-user li120415323@gmail.com

# Load model
module load bioinfo-tools BEDTools/2.27.1
BED_PATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllBED
#numSample=5
# bin option 1:
WindowSize=30
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}
echo "Initiate and clean the folder..."
mkdir -p $OutPATH
rm -f $OutPATH/*
echo "Getting all the records into one BED file..."
ls ${BED_PATH}/*General.bed > ${OutPATH}/Samples.list # get the sample names
cat ${BED_PATH}/*General.bed > ${OutPATH}/All.tmp.bed
sort -k1,1 -k2,2n ${OutPATH}/All.tmp.bed | bedtools merge -i stdin > ${OutPATH}/All.tmp.sort.merge.bed
# bin option 1 main:
echo "Generating window fo size $WindowSize across the genome..."
bedtools makewindows -b  ${OutPATH}/All.tmp.sort.merge.bed -w $WindowSize | awk '{OFS="\t"; print $0, $3-$2}' > ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
echo "Finished, the number of windows in the file is:"
wc -l ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
echo "Bin all Samples, WindowSize: $WindowSize"
#cat ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.bed | awk '{OFS="\t"; print $0, $3-$2}'  > ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
cat ${OutPATH}/Samples.list | \
while read p; do
  BaseName=$(basename ${p%%.*})
  echo "Sample name: $BaseName, Size:"
  bedtools intersect -wao -f 0.5 -a ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed -b $p > ${OutPATH}/${BaseName}.binWin${WindowSize}.bed
  wc -l ${OutPATH}/${BaseName}.binWin${WindowSize}.bed
done


#sbatch /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/4_binVariant_bin30bp.sh > /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/4_binVariant_bin30bp.log
#sbatch /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/4_binVariant_bin10bp.sh > /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/4_binVariant_bin10bp.log

# bin option 2:
WindowSize=10
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}
echo "Initiate and clean the folder..."
mkdir -p $OutPATH
rm -f $OutPATH/*
echo "Getting all the records into one BED file..."
ls ${BED_PATH}/*General.bed > ${OutPATH}/Samples.list # get the sample names
cat ${BED_PATH}/*General.bed > ${OutPATH}/All.tmp.bed
sort -k1,1 -k2,2n ${OutPATH}/All.tmp.bed | bedtools merge -i stdin > ${OutPATH}/All.tmp.sort.merge.bed
#NOTE: wc -l /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/All.tmp.sort.merge.bed #=> 30598
# bin option 2, main:
echo "Generating window fo size $WindowSize across the genome..."
bedtools makewindows -b  ${OutPATH}/All.tmp.sort.merge.bed -w $WindowSize | awk '{OFS="\t"; print $0, $3-$2}' > ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
echo "Finished, the number of windows in the file is:"
wc -l ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
echo "Bin all Samples, WindowSize: $WindowSize"
#cat ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.bed | awk '{OFS="\t"; print $0, $3-$2}'  > ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
cat ${OutPATH}/Samples.list | \
while read p; do
  BaseName=$(basename ${p%%.*})
  echo "Sample name: $BaseName, Size:"
  bedtools intersect -wao -f 0.5 -a ${OutPATH}/All.tmp.sort.merge_win${WindowSize}.size.bed -b $p > ${OutPATH}/${BaseName}.binWin${WindowSize}.bed
  wc -l ${OutPATH}/${BaseName}.binWin${WindowSize}.bed
done

#cound NA in the population:
#cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10
#grep 'NA' All.tmp.bed  | wc -l
#85397
#wc -l All.tmp.bed
#4602489

# Compare the line in the file, need to dedup
#wc -l /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/All.tmp.sort.merge_win10.size.bed
#22235594 /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/All.tmp.sort.merge_win10.size.bed

#just testing 50 window table size: cat ${BED_PATH}/*sort.rmOverlap.bed > ${OutPATH}/All.tmp.bed
#bedtools makewindows -b  BinVariants/All.tmp.sort.merge.bed -w 50 > BinVariants/All.tmp.sort.merge_win50.bed # old scripts
# size: 2899352 BinVariants/All.tmp.sort.merge_win50.bed; CNVnator raw window number: 2166502

#test:
# TEST LOG! aim to have same funtion as CNVnator's script: 3_CNV_matrix.sh:
# Debuging
# cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants
# # get all
# cat  testHead2.bed | bedtools makewindows -b stdin -w 50 | awk '{print $0,$3-$2}' | less -S
# sbatch -A sens2016007 -p core -n 8 -t 4:00:00 --mail-user li120415323@gmail.com -J testDifferentBins
