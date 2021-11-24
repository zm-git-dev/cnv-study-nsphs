#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2-15:00:00
#SBATCH -J ${list_1##*/}_CNVnatorPartition
echo "Starting at:"
date
# load the ROOT library
module load ROOT/6.06.08
module load bioinfo-tools CNVnator
CNVnatorFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_${list_1##*/}
Inter_OutputFolder=\${CNVnatorFolder}/IntermediaOutput
Log_Folder=\${CNVnatorFolder}/Log
CNVnator_1_BinSizeList=\${CNVnatorFolder}/1_BinSizeList.txt
# define the SN for analysis, to reduce the warning messages for e.g SN: 'GL000223.1'
chrom="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
while read RootPath bin_size
do
  BaseName=\${RootPath##*/}
  LogOutput=\${Log_Folder}/\${BaseName%%.*}.log
  echo Analysing \${BaseName}
  date
  # >>>RD SIGNAL PARTITIONING
  cnvnator -root \$RootPath -partition \$bin_size -chrom \$chrom >> \$LogOutput
  echo Finish \${BaseName} at:
  date