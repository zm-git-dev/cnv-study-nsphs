#!/bin/bash -l

#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:45:00
#SBATCH -J PasteBin30bp
#SBATCH --mail-type=ALL
#SBATCH --mail-user li120415323@gmail.com

WindowSize=30
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
for list in lists*; do echo "Processing $list"; paste $(awk '{OFS="\t"; print $12}' $list) > merge${list##lists}; done
paste merge* > Paste_${WindowSize}.txt

# sbatch /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/5_PasteBin30.sh > /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/5_PasteBin30.log
# ! if start from the same folder, slurm.out file will be removed!

#for list in head10_lists10; do echo "Processing $list"; paste $(cat $list) -d "," > Test_merge; done
# for Bin 30 bp:
# Check if all the files have the same number of row:
# cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10
# ls ./*.binWin10.bed | wc -l
#
# module load R/3.4.3
# module load R_packages/3.4.3
