module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3

## Format CNV matrix for analysis
# test loading data with data table
# rm(list=ls())
# setwd("/castor/project/proj_nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder")
# library(data.table)
# my_file <- list.files(pattern = "Population_collapse*", full.names=TRUE)
# my_table <- fread(file = my_file)
# #change . to 2 copies, check if NA is still here:
# summary(my_table==".")
# my_table[my_table == "."] <- 2
# summary(my_table==".")
# # save R data
# saveRDS(my_table,"Format_populationMatrix_bin10.rds")
# # new environment:
# CNVbed <- readRDS("/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Format_populationMatrix_bin10.rds")
# rownames(CNVbed) <- paste0(CNVbed$Chr, ":", CNVbed$Start, "-", CNVbed$End)
# saveRDS(CNVbed,"Format_population_bin10_rowname.rds")

# in BASH:

cd /home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/7_subScripts
for i in {1..22}
do
   cat > CNV_Pea3_Chr${i}.sh <<- EOM
#!/bin/bash -l
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 9:00:00
#SBATCH -J CNV_Pea3_Chr${i}
#SBATCH --mail-type=ALL
#SBATCH --mail-user li120415323@gmail.com

module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
R CMD BATCH --no-restore --no-save CNV_Pea3_Chr${i}.R
EOM
done

for i in {1..22}
do
   cat > CNV_Pea3_Chr${i}.R <<- EOM
source("/home/zhiwei94/CNVscripts/Manta1_5Folder/FormatPopulationScripts/7_WholeG_LinearRegression.R")
LR_function(pea_3, CNVbed, AgeSex, ${i})
EOM
done

for each in *.sh; do echo "sbatch ${each} > ${each%%.sh}.log"; done
