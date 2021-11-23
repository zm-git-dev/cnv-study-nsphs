# List of input files:
## List of BAM files, including the soft link, saved to home direcotry
ll /proj/sens2016007/nobackup/bam-files/*.bam > ~/ListOfBamFiles.txt
cat ~/ListOfBamFiles.txt | awk '{print $9}' > ~/NoSoftLinkListOfBamFiles.txt
## Split the 1047 sample list into 100 files, ~/NoSoftLinkListOfBamFiles.txt
### Go to the directory for saving list
cd ~/BamList
## Group the 1047 bam files by 100 sample size:
split -l100 ~/NoSoftLinkListOfBamFiles.txt --verbose
## Save the list files into one single list: ListOfGroup
ll ~/BamList/x* | awk '{print $9}' > ~/BamList/ListOfGroup.txt

# CNVnator:1:Generating root files:
## cat /home/zhiwei94/CNVscripts/CNVnatorFolder/RootMain/root_all_generate.sh
### Script saved to:
cd /home/zhiwei94/CNVscripts/CNVnatorFolder/RootMain/AllScript
### Use loop to generat 11 script for Root test
### Inside the loop function: ~/BamList/ListOfGroup.txt:path of the 11 group, $list_1:each group
while read list_1
do
  cat > ${list_1##*/}_root.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 3-12:00:00
#SBATCH -J ${list_1##*/}_CNVnatorRoot
echo "Starting at:"
date
# load the ROOT library
module load ROOT/6.06.08
module load bioinfo-tools CNVnator
CNVnatorFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_${list_1##*/}
Inter_OutputFolder=\${CNVnatorFolder}/IntermediaOutput
Log_Folder=\${CNVnatorFolder}/Log
# Create the output folder, if not exist
mkdir -p \$Inter_OutputFolder
mkdir -p \$Log_Folder
# define the SN for analysis, to reduce the warning messages for e.g SN: 'GL000223.1'
chrom="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
EOM
  while read BamFile
  do
    BaseName=${BamFile##*/}
    cat >> ${list_1##*/}_root.sh <<- EOM
echo Reading BAM file ${BaseName%%.*}
cd \$Inter_OutputFolder
RootOutput=\${Inter_OutputFolder}/${BaseName/bam/root}
echo Extracting reads...
cnvnator -root \$RootOutput -chrom \$chrom -tree $BamFile -unique > \${Log_Folder}/${BaseName%%.*}.log
EOM
  done < $list_1
done < ~/BamList/ListOfGroup.txt
###########################################
## Run all the root job
for each in *; do echo "sbatch ${each} > ${each%.*}.log"; done
sbatch xaa_root.sh > xaa_root.log
sbatch xab_root.sh > xab_root.log
sbatch xac_root.sh > xac_root.log
sbatch xad_root.sh > xad_root.log
sbatch xae_root.sh > xae_root.log
sbatch xaf_root.sh > xaf_root.log
sbatch xag_root.sh > xag_root.log
sbatch xah_root.sh > xah_root.log
sbatch xai_root.sh > xai_root.log
sbatch xaj_root.sh > xaj_root.log
sbatch xak_root.sh > xak_root.log
###########################################

# CNVnator:2:Computing optimal binsize for each sample:
## Script to generate all binSize scripts for xaa~xak
## cat /home/zhiwei94/CNVscripts/CNVnatorFolder/MainBinSize/binSize_all_generate.sh
### Script saved to:
cd /home/zhiwei94/CNVscripts/CNVnatorFolder/MainBinSize/AllScript
# Use loop to generat 11 script for BinSize processing
# General info:
while read list_1
do
  cat > ${list_1##*/}_BinSize.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-08:00:00
#SBATCH -J ${list_1##*/}_CNVnatorBinSize
echo "Starting at:"
date
# load the ROOT library
module load ROOT/6.06.08
module load bioinfo-tools CNVnator
declare -a bin_size_List=(70 85 100 150 200 250)
CNVnatorFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_${list_1##*/}
Inter_OutputFolder=\${CNVnatorFolder}/IntermediaOutput
Log_Folder=\${CNVnatorFolder}/Log
# Create the output folder, if not exist
mkdir -p \$Inter_OutputFolder
mkdir -p \$Log_Folder
# The results of binSize testing will be the List1: optimal binSize for each BAM and the List2: BAMs that did not find optimal binSize in the given list
OutputList=\${CNVnatorFolder}/1_BinSizeList.txt
AbnormalList=\${CNVnatorFolder}/2_AbnormalList.txt
# Remove the previous binSize testing output
rm -f \$OutputList
rm -f \$AbnormalList
# Define the SN for analysis, to reduce the warning messages for e.g SN: 'GL000223.1'
chrom="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
while read BamFile
do
  BaseName=\${BamFile##*/}
  echo Reading BAM file \${BaseName%%.*}
  cd \$Inter_OutputFolder
  RootOutput=\${Inter_OutputFolder}/\${BaseName/bam/root}
  tmpOutput=\${Inter_OutputFolder}/\${BaseName/bam/tmp}.binStat.txt
  #remove previous tmpOutput
  rm -f \$tmpOutput
  for bin_size in "\${bin_size_List[@]}"
  do
    echo Calculating stat for bin size \$bin_size >> \${Log_Folder}/\${BaseName%%.*}.log
    # >>GENERATING A HISTOGRAM
    echo GENERATING A HISTOGRAM >> \${Log_Folder}/\${BaseName%%.*}.log
    # only specify the folder of all single reference sequences, not each.
    cnvnator -root \$RootOutput -his \$bin_size -chrom \$chrom -d /proj/sens2016007/nobackup/Reference/ >> \${Log_Folder}/\${BaseName%%.*}.log
    # >>>CALCULATING STATISTICS
    echo CALCULATING STATISTICS >> \${Log_Folder}/\${BaseName%%.*}.log
    # For the stat output, if the mean\$7/std\$9 ratio is outside (4~5) from autosome, add "False" flag in column 4, otherwise add "True" flag
    cnvnator -root \$RootOutput -stat \$bin_size -chrom \$chrom | tail -n2 | head -n1 | awk -v my_var1=\${BaseName/bam/root} -v my_var2=\$bin_size -v OFS="\t" '{if (\$7/\$9>4&&\$7/\$9<5) print my_var1, my_var2, \$7/\$9, "True"; else print my_var1, my_var2, \$7/\$9, "False"}' | sort -k3n,3 >> \$tmpOutput
  done
  # Reverse sort the flag column, true will popup to the first row
  if [ "\$(sort -rk4,4 \$tmpOutput | awk 'NR==1 {print \$4}')" == "True" ]
  then
    # Save the .root files full PATH, selected bin_size and Average RD's mean/stdev ratio closest to 4 for downstream analyses
    echo \${BaseName%%.*} has bin_size showing good mean/stdev ratio >> \${Log_Folder}/\${BaseName%%.*}.log
    awk -v my_var1=\$RootOutput -v OFS="\t" '{if (\$4=="True") print my_var1, \$2, \$3, \$4}' \$tmpOutput | sort -k3n,3 | sort -k1,1 -u >> \$OutputList
  else
    #echo \$(sort -rk4,4 \$tmpOutput | awk 'NR==1 {print \$4}')
    echo No suitable bin_size selected for file \${BaseName%%.*} >> \${Log_Folder}/\${BaseName%%.*}.log
    echo STATISTICS results saved to file: \$AbnormalList >> \${Log_Folder}/\${BaseName%%.*}.log
    # Save the PATH to the .BAM files for re-calculating the insert size STATISTICS
    awk -v File_Path=\$BamFile -v OFS="\t" '{print File_Path, \$2, \$3, \$4}' \$tmpOutput >> \$AbnormalList
  fi
done < $list_1
EOM
done < ~/BamList/ListOfGroup.txt
###########################################
# Submitting all the jobs
for each in x*; do echo "sbatch ${each} > ${each%.*}.log"; done
###########################################
# There are 4 samples from 2 subsets did not get good bin size, generate two jobs for the two subset to rerun with more binSize
# Get the sample names by:
InputBamList=$(awk -v OFS="\t" '{print $1}' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xah/2_AbnormalList.txt | sort -u)
InputBamList=$(awk -v OFS="\t" '{print $1}' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xaj/2_AbnormalList.txt | sort -u)
# Testing loop:
for BamFile in $InputBamList; do echo $BamFile; done
# rerun the jobs, there were 4 samples with no good bin size
sbatch xah_BadBin_rerun.sh > xah_BadBin.log
sbatch xaj_BadBin_rerun.sh > xaj_BadBin.log
# The re-run scripts is in the MainBinSize/AllScript folder as: xah_CNVnator_binSize and xaj_CNVnator_binSize
# Mean value of optimal bin_size of all samples:
awk '{sum+= $2}END{print "Average: " sum/NR }' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/1_BinSizeList.txt
#Average: 91.533
###########################################
# Check if the finished job contain error messages
# search Ignore case distinctions
grep -i 'ERROR' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xak/Log/*
cat /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/Log/* | grep -i  'ERROR'
ll /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/Log | grep 'total'
ls -d /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/Log/* | less -S
###########################################

# CNVnator:3:PARTITIONING
## The read depth signals are 'smoothen' by partitioning algorithm,
## For this step, only the list of the group is needed, it get the .root bin_size form 1_BinSizeList.txt of each CNVnator folder
## Generate PARTITIONING job scripts for ALL BAM subset
## cat /home/zhiwei94/CNVscripts/CNVnatorFolder/MainPartition/Partition_all_generate.sh
### Script saved to:
SavePATH=/home/zhiwei94/CNVscripts/CNVnatorFolder/MainPartition/AllScript
mkdir -p $SavePATH
cd $SavePATH
while read list_1
do
  cat > ${list_1##*/}_Partition.sh <<- EOM
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
done < \$CNVnator_1_BinSizeList
EOM
done < ~/BamList/ListOfGroup.txt
###########################################
## Run all the root job
for each in *; do echo "sbatch ${each} > ${each%.*}.log"; done
sbatch xaa_Partition.sh > xaa_Partition.log
sbatch xab_Partition.sh > xab_Partition.log
sbatch xac_Partition.sh > xac_Partition.log
sbatch xad_Partition.sh > xad_Partition.log
sbatch xae_Partition.sh > xae_Partition.log
sbatch xaf_Partition.sh > xaf_Partition.log
sbatch xag_Partition.sh > xag_Partition.log
sbatch xah_Partition.sh > xah_Partition.log
sbatch xai_Partition.sh > xai_Partition.log
sbatch xaj_Partition.sh > xaj_Partition.log
sbatch xak_Partition.sh > xak_Partition.log
###########################################

# CNVnator:4:CNVs calling
## For this step, only the list of the group is needed, it get the .root bin_size form 1_BinSizeList.txt of each CNVnator folder
## Generate PARTITIONING job scripts for ALL BAM subset
## cat /home/zhiwei94/CNVscripts/CNVnatorFolder/MainCallCNVs/CallCNVs_all_generate.sh
### Script saved to:
SavePATH=/home/zhiwei94/CNVscripts/CNVnatorFolder/MainCallCNVs/AllScript
mkdir -p $SavePATH
cd $SavePATH
# General info:
while read list_1
do
  cat > ${list_1##*/}_CallCNVs.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:20:00
#SBATCH -J ${list_1##*/}_CNVnatorCallCNVs
echo "Starting at:"
date
# load the ROOT library
module load ROOT/6.06.08
module load bioinfo-tools CNVnator
CNVnatorFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_${list_1##*/}
Inter_OutputFolder=\${CNVnatorFolder}/IntermediaOutput
Log_Folder=\${CNVnatorFolder}/Log
VariantPath=\${CNVnatorFolder}/3_Variant_output
mkdir -p \$VariantPath
CNVnator_1_BinSizeList=\${CNVnatorFolder}/1_BinSizeList.txt
# define the SN for analysis, to reduce the warning messages for e.g SN: 'GL000223.1'
chrom="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
rm -f \$AllLog
while read RootPath bin_size
do
  BaseName=\${RootPath##*/}
  LogOutput=\${Log_Folder}/\${BaseName%%.*}.log
  echo Call CNVs for \${BaseName} >> \$LogOutput
  SECONDS=0
  # >>>CNV CALLING, calls are printed to STDOUT. TODO: save the VCF.txt to output home directory or a seperate VCF file?
  cnvnator -root \$RootPath -call \$bin_size > \${VariantPath}/\${BaseName/root/variant.txt}
  echo Finish \${BaseName} at:
  date
  duration=\$SECONDS
  echo "For \${BaseName%%.*} CNVsCalling, \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."
done < \$CNVnator_1_BinSizeList
EOM
done < ~/BamList/ListOfGroup.txt
###########################################
for each in *; do echo "sbatch ${each} > ${each%.*}.log"; done
###########################################

# CNVnator:5:Formating stdout CNV results: see CNV_matrix.sh for more detail
## All txt results
OutPATH=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSizeCV
mkdir -p $OutPATH
VariantTXT=$(ls -d  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/3_Variant_output/*)
for CNVnatorCall in ${VariantTXT}
do
  BaseName=${CNVnatorCall##*/}
  awk '{if ($9>=0&&$9<=0.5&&$3>=1000&&$5<=0.00005) print $0}' $CNVnatorCall | awk -v OFS="\t" '{split($2, a, /:|-/); print a[1], a[2]-1, a[3], $4*2, $5, $6, $9, $3}' \
  > ${OutPATH}/${BaseName%.txt}_selected.bed
done
# also filter e-val2? what is e-val2?
awk '{if ($6<=0.00005) print $0}'  /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSizeCV/P4604_170.clean.dedup.variant_selected.bed | less -S
awk '{if ($6<=0.05&&$4<=1) print $0}'  /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSizeCV/P4604_170.clean.dedup.variant_selected.bed | wc -l
###########################################


###########################################
# Manta 1.5.0
## There is mainly one step in Manta variants calling step, reporting SVs including deletions, duplications, small indels, inversions
###########################################
SavePATH=/home/zhiwei94/CNVscripts/Manta1_5Folder/AllScriptCores
mkdir -p $SavePATH
cd $SavePATH
# Loop for creating Manta variants calling scripts for each group
while read list_1
do
  cat > ${list_1##*/}_MantaCores.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 6-18:00:00
#SBATCH -J ${list_1##*/}_Manta1.5_8core
echo "Starting at:"
date
# load MANTA
#module load bioinfo-tools manta/1.0.3
# Local MANTA1.5
MantaFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_${list_1##*/}
ReferenceGenome=/proj/sens2016007/nobackup/Reference/human_g1k_v37.fasta
MANTA_INSTALL_PATH=/proj/sens2016007/nobackup/Tools/Manta_1.5.0/
Inter_OutputFolder=\${MantaFolder}/IntermediaOutput
LogFolder=\${MantaFolder}/Logs
mkdir -p \$LogFolder
AllLog=\${MantaFolder}/All.log
rm -f \$AllLog
while read BamFile
do
  BaseName=\${BamFile##*/}
  MANTA_ANALYSIS_PATH=\${MantaFolder}/IntermediaOutput/\${BaseName%%.*}_Folder
  LogOutput=\${LogFolder}/\${BaseName%%.*}.log
  rm -rf \$MANTA_ANALYSIS_PATH
  mkdir -p \$MANTA_ANALYSIS_PATH
  echo Analysing \${BaseName}
  echo Analysing \${BaseName} >> \$AllLog
  date >> \$AllLog
  date
  SECONDS=0
  cd \$MANTA_ANALYSIS_PATH
  \${MANTA_INSTALL_PATH}/bin/configManta.py --bam \$BamFile --referenceFasta \$ReferenceGenome --runDir \${MANTA_ANALYSIS_PATH} > \$LogOutput
  \${MANTA_ANALYSIS_PATH}/runWorkflow.py --mode local -j 8 >> \$LogOutput
  echo done >> \$AllLog
  date >> \$LogOutput
  duration=\$SECONDS
  echo "For \${BaseName%%.*}, \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed." >> \$AllLog
  echo Finish analysis for \${BaseName}
done < $list_1
EOM
done < ~/BamList/ListOfGroup.txt
###########################################
# submit jobs, try xaa first.
## xak subset runtime change to 4-18:00:00
for each in xa*.sh; do echo "sbatch  $each > ${each%.sh}.id "; done
## submitted!
###########################################
# Check error:
## One file no result:
ls /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaj/IntermediaOutput/*Folder/results/variants/diploidSV.vcf.gz | wc -l
##99
grep -i 'error' /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaj/IntermediaOutput/*_Folder/workflow.error.log.txt
## File did not complete: P4605_224
less -S /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaj/IntermediaOutput/P4605_224_Folder/workflow.error.log.txt
grep 'P4605_224' ~/NoSoftLinkListOfBamFiles.txt > /home/zhiwei94/BamList/Manta1_5_rerun
ls -d /home/zhiwei94/BamList/Manta1_5_rerun > ~/BamList/ListManta1_5.txt
## Rerun P4605_224
SavePATH=/home/zhiwei94/CNVscripts/Manta1_5Folder/AllScriptCores
mkdir -p $SavePATH
cd $SavePATH
while read list_1
do
  cat > ${list_1##*/}_MantaCores.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J ${list_1##*/}_Manta1.5_8core
echo "Starting at:"
date
# load MANTA
#module load bioinfo-tools manta/1.0.3
# Local MANTA1.5
MantaFolder=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_${list_1##*/}
ReferenceGenome=/proj/sens2016007/nobackup/Reference/human_g1k_v37.fasta
MANTA_INSTALL_PATH=/proj/sens2016007/nobackup/Tools/Manta_1.5.0/
Inter_OutputFolder=\${MantaFolder}/IntermediaOutput
LogFolder=\${MantaFolder}/Logs
mkdir -p \$LogFolder
AllLog=\${MantaFolder}/All.log
rm -f \$AllLog
while read BamFile
do
  BaseName=\${BamFile##*/}
  MANTA_ANALYSIS_PATH=\${MantaFolder}/IntermediaOutput/\${BaseName%%.*}_Folder
  LogOutput=\${LogFolder}/\${BaseName%%.*}.log
  rm -rf \$MANTA_ANALYSIS_PATH
  mkdir -p \$MANTA_ANALYSIS_PATH
  echo Analysing \${BaseName}
  echo Analysing \${BaseName} >> \$AllLog
  date >> \$AllLog
  date
  SECONDS=0
  cd \$MANTA_ANALYSIS_PATH
  \${MANTA_INSTALL_PATH}/bin/configManta.py --bam \$BamFile --referenceFasta \$ReferenceGenome --runDir \${MANTA_ANALYSIS_PATH} > \$LogOutput
  \${MANTA_ANALYSIS_PATH}/runWorkflow.py --mode local -j 8 >> \$LogOutput
  echo done >> \$AllLog
  date >> \$LogOutput
  duration=\$SECONDS
  echo "For \${BaseName%%.*}, \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed." >> \$AllLog
  echo Finish analysis for \${BaseName}
done < $list_1
EOM
done < ~/BamList/ListManta1_5.txt
