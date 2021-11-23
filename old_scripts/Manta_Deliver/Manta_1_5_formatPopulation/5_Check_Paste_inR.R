# for bin=10
# cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10
# ls ./*.binWin10.bed | wc -l
#
# module load R/3.4.3
# module load R_packages/3.4.3
library(data.table)
my_file_list <- list.files(pattern = "*.binWin10.bed", full.names=TRUE)
# dt <- fread(my_file_list)
# dt <- fread("ls ./*.binWin10.bed")
# took a long time:
# use first 5 file for testing
#https://stackoverflow.com/questions/42819343/adding-a-column-characters-based-on-file-name-in-r
# https://stackoverflow.com/questions/45514204/how-to-add-a-index-by-set-of-data-when-using-rbindlist/45515611#45515611
my_table <- rbindlist(lapply(my_file_list[1:2], fread, header = FALSE, select = c(1:4,12)), idcol = "FileName")
#rename_my_table <- my_table[, FileName := factor(FileName, labels = basename(my_file_list[1:2])) ]
rename_my_table <- my_table[, FileName := basename((my_file_list)[FileName]) ]
