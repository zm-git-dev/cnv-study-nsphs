#!/bin/bash
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllVCF
mkdir -p $OutPATH
# save the VCF files directory to one txt file.
ls -d  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_*/IntermediaOutput/*_Folder/results/variants/diploidSV.vcf.gz > $(dirname $OutPATH)/sample_files_VCF.txt
#VariantGZ=$(ls -d  /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_*/IntermediaOutput/*_Folder/results/variants/diploidSV.vcf.gz)
# for Call in ${VariantGZ}
# do
#   BaseName=$(basename $(echo ${Call} | sed 's/.*\(IntermediaOutput\)/\1/g' | cut -d'/' -f-2) )
#   echo Saving to ${OutPATH}/${BaseName%_Folder}.PASS.bed
#   zcat $Call | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe | awk '(($11=="DUP"||$11=="DEL")&&$12=="PASS"&&$1==$4&&$6>$2){OFS="\t";print $1,$2,$6,$6-$2,$0}(($11=="DUP"||$11=="DEL")&&$12=="PASS"&&$1==$4&&$6<$2){OFS="\t";print $1,$5,$3,$3-$5,$0}' > ${OutPATH}/${BaseName%_Folder}.PASS.bed
# done
#loop through lines in the VCF.txt fiile
while read p; do
  BaseName=$(basename $(echo ${p} | sed 's/.*\(IntermediaOutput\)/\1/g' | cut -d'/' -f-2) )
  echo copy "$p" to $OutPATH
  cp $p ${OutPATH}/${BaseName%_Folder}.vcf.gz
done <$(dirname $OutPATH)/sample_files_VCF.txt

# remove the 10 samples mentioned in /proj/sens2016007/nobackup/work/analyis/ASA_joint_calling_2017_04_03/exlcluded_samples_do_to_genotype_conc.txt
while read p; do
  echo remove file: ${p}.vcf.gz
  rm ${OutPATH}/${p}.vcf.gz
done < /proj/sens2016007/nobackup/work/analyis/ASA_joint_calling_2017_04_03/exlcluded_samples_do_to_genotype_conc.txt

# Remove list: Collected in /proj/sens2016007/nobackup/diana/data/excludeSamples.txt
while read p; do
  echo remove file: ${p}.vcf.gz
  rm ${OutPATH}/${p}.vcf.gz
done < /proj/sens2016007/nobackup/diana/data/excludeSamples.txt

#check the number of vcf.gz files:
echo Finish. Number of samples in the folder: $(ls $OutPATH | wc -l)
