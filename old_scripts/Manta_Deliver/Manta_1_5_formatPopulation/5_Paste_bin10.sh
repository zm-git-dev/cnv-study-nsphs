#!/bin/bash -l

#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:45:00
#SBATCH -J PasteBin10bp
#SBATCH --mail-type=ALL
#SBATCH --mail-user li120415323@gmail.com

WindowSize=10
BED_PATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}/PasteFolder
# bedtools merge -a file:
BEDtools_a=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}/All.tmp.sort.merge_win${WindowSize}.size.bed
# Check if each BED has the same rows as BEDtools_a file. Generate files(lists) with only CN estimation
echo "Checking number of rows and generate list format CN info"
for BED_file in ${BED_PATH}/*.binWin${WindowSize}.bed
  do if [ "$(wc -l < $BED_file)" -eq "$(wc -l < $BEDtools_a)" ]
  then
    #echo "$BED_file Pass"
    awk '{print $12}' $BED_file > ${BED_file%.bed}.list
  else
    echo "File $BED_file don't have the same rows as BEDtools_a file, exiting..."
    exit
  fi
done
echo "Initiate and clean the folder..."
mkdir -p $OutPATH
rm -f $OutPATH/*
# go to the PasteFolder:
cd $OutPATH
ls -1 ${BED_PATH}/*.binWin${WindowSize}.list > ListFile.txt
ls -1 ${BED_PATH}/*.binWin${WindowSize}.list | split -l 100 -d - lists
# use paste with awk to speed up the pasting process e.g: paste *.txt > combined.txt
for list in lists*; do echo "Processing $list"; paste $(awk '{OFS="\t"; print $12}' $list) > merge${list##lists}; done
paste merge* > Paste_${WindowSize}.txt

# the same list loop order to generate the sample name list:
cd $OutPATH
suffix=".binWin${WindowSize}.list"
# remove file path and suffix:
for list in lists*; do cat ${list} | sed -e 's!.*/!!' -e "s/$suffix$//"  ; done > Sample_name.txt
# change name: would it be more efficient using dictionary function in Python?
while read each id others;
do
  BaseName=${each##*/}
  echo -e "${BaseName%%.*}\t$id"
done < /proj/sens2016007/nobackup/NSPHS_phenotype_data/WGS_kodnyckel > namesDic.txt
# Convert to other name:
for each in $(cat Sample_name.txt);
do
  #echo "Filename $each"
  while read Filename Id;
  do
    if [ "$each" = "$Filename" ]
    then
      echo "$Id"
    fi
  done < namesDic.txt
done > Sample_name_changeName.txt
#add the name to the matrix
paste -sd"\t" Sample_name_changeName.txt | cat - Paste_${WindowSize}.txt > Paste_${WindowSize}.withHeader.bed
# check if the number of columns is fine?
# awk '{print NF}' Paste_${WindowSize}.withHeader.bed | uniq -c
# 22235595 1021
# get the BED coordiante info: Chr, Start, End. Not including size?
WindowList=${BED_PATH}/All.tmp.sort.merge_win${WindowSize}.size.bed
echo -e "Chr\tStart\tEnd" > WindowList.txt
awk '{OFS="\t"; print $1,$2,$3}' $WindowList >> WindowList.txt
# Paste population with the window list
paste WindowList.txt Paste_${WindowSize}.withHeader.bed > Population_${WindowSize}.withHeader.bed
# collpase the populstion matrix:

###test the hOSCAR regioin reported in CNVnator:
#19	54555500	54560700
#echo -e "19\t54555500\t54560700\t" > hOSCAR_coor.bed # remove the last tab!
# us vim to remove the last tab in the hOSCAR file
#module load bioinfo-tools BEDTools/2.27.1
#bedtools intersect -wao -a hOSCAR_coor.bed -b Population_${WindowSize}.withHeader.bed | less -S
###
#
# sbatch /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/5_PasteBin10.sh > /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/5_PasteBin10.log
# ! if start from the same folder, slurm.out file will be removed!

#for list in head10_lists10; do echo "Processing $list"; paste $(cat $list) -d "," > Test_merge; done
# for Bin 30 bp:
# Check if all the files have the same number of row:
# cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10
# ls ./*.binWin10.bed | wc -l
#
# module load R/3.4.3
# module load R_packages/3.4.3
