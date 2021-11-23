#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 3:00:00
#SBATCH -J collapseMatrix
echo "Starting at:"
date
WindowSize=10
#BED_PATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize${WindowSize}/PasteFolder
cd $OutPATH
# define input population matrix formatted by window size:
#inBed="chr19_54399999_54599999.bed"
# test with first 200 lines:
head -n200 Population_${WindowSize}.withHeader.bed > head200_Population_${WindowSize}.withHeader.bed
#inBed=Population_${WindowSize}.withHeader.bed
inBed=head200_Population_${WindowSize}.withHeader.bed
Lines=$(wc -l < $inBed) #! $echo "$Lines" and $echo $Lines give different values
# The population matirx in read into 4 field: Chromosome, StartPosition, EndPosition, Population_genotype(1021 genotype, read as a string)
PreChrom=1
PreStart=0
PreEnd=0
PreGenotype=""
counter=0
while read Chrom Start End Genotype;
do
  counter=$((counter+1))
  if [ "$counter" = $Lines ] # If read to the last line, [ "$counter" = "$Lines "] doesn't work, Lines=$(wc -l < $inBed) causes error when "$Lines"
  then
    #echo "Last line"
    if [ "$Chrom" = "$PreChrom" ] && [ "$Start" = "$PreEnd" ] && [ "$Genotype" = "$PreGenotype" ]
    then
      PreEnd="$End"
      echo -e "$PreChrom\t$PreStart\t$PreEnd\t$PreGenotype"
    else
      echo -e "$PreChrom\t$PreStart\t$PreEnd\t$PreGenotype"
      echo -e "$Chrom\t$Start\t$End\t$Genotype"
    fi
  elif  [ "$Chrom" = "$PreChrom" ] && [ "$Start" = "$PreEnd" ] && [ "$Genotype" = "$PreGenotype" ] # if the current window has the same genotype as the previous window:
  then
    PreEnd="$End"
  else
    if [ "$counter" != 1 ]
    then
     echo -e "$PreChrom\t$PreStart\t$PreEnd\t$PreGenotype"
    fi
    PreChrom="$Chrom"
    PreStart="$Start"
    PreEnd="$End"
    PreGenotype="$Genotype"
  fi
done < $inBed > Population_collapse_bin${WindowSize}.bed
## Result:
wc -l $inBed
##Manta:22235595 CNVnator:2166502
wc -l Population_collapse_bin${WindowSize}.bed
##Manta:124527 CNVnator:263831

module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3

## Format CNV matrix for analysis in R
#test loading data with data table
rm(list=ls())
setwd("/castor/project/proj_nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder")
library(data.table)
my_file <- list.files(pattern = "Population_collapse*", full.names=TRUE)
my_table <- fread(file = my_file)
#change . to 2 copies, check if NA is still here:
summary(my_table==".")
my_table[my_table == "."] <- 2
summary(my_table==".")
# save R data
saveRDS(my_table,"Format_populationMatrix_bin10.rds")
# new environment:
CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_populationMatrix_bin10.rds")
rownames(CNVbed) <- paste0(CNVbed$Chr, ":", CNVbed$Start, "-", CNVbed$End)
saveRDS(CNVbed,"Format_population_bin10_rowname.rds")
