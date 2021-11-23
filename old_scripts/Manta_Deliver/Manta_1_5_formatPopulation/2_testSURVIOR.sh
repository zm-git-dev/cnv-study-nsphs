#!/bin/bash -l

#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:45:00
#SBATCH -J SURVIVOR_PopVCF
#SBATCH --mail-type=ALL
#SBATCH --mail-user li120415323@gmail.com

module load bioinfo-tools SURVIVOR/1.0.3
module load bioinfo-tools picard/2.20.4
#module load bioinfo-tools GATK/4.1.1.0
module load bioinfo-tools bcftools/1.8
#/sw/bioinfo/GATK/4.1.1.0/bianca/gatk
OutPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/SURVIRO/samples
PopulationPath=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/SURVIRO
VcfPATH=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllVCF
PICARD="/sw/bioinfo/picard/2.20.4/bianca/picard.jar"
mkdir -p $OutPATH
# clean the files in the directory.
rm -f ${OutPATH}/*
# folder with vcf.gz files: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/AllVCF
while read each id others;
do
  BaseName=${each##*/}
  echo -e "${BaseName%%.*}\t$id"
done < /proj/sens2016007/nobackup/NSPHS_phenotype_data/WGS_kodnyckel > ${OutPATH}/namesDic.txt
# List all samples, and unzip them:
ls -d $VcfPATH/* > ${OutPATH}/AllSamples.txt
# read all the files, unzip
while read p; do
  #TODO: filter variants and change name with RenameSampleInVcf (Picard)??
  BaseName=$(basename $p)
  #echo "$BaseName"
  # change sample names:
  while read Filename Id;
  do
    if [ "${BaseName%%.*}" = "$Filename" ]
    then
      echo -e "$BaseName\t$Id"
      java -jar $PICARD RenameSampleInVcf \
      INPUT=$p \
      OUTPUT=${OutPATH}/${BaseName%%.*}_${Id}.vcf \
      NEW_SAMPLE_NAME=$Id
      bcftools view -f PASS ${OutPATH}/${BaseName%%.*}_${Id}.vcf > ${OutPATH}/${BaseName%%.*}_${Id}_PASS.vcf
    fi
  done < ${OutPATH}/namesDic.txt
  # keep only PASS
  #zcat $p > ${BaseName%.vcf.gz}.vcf
  #echo $BaseName
done < ${OutPATH}/top5.txt
#test all PASS: awk -F '\t' '!($0 ~ /^\#/){print $7}' KA09-155_KA09-0155_PASS.vcf | sort | uniq -c
ls -d ${OutPATH}/*PASS.vcf > ${OutPATH}/sample_files.txt

# test running SURVIVOR
#module load bioinfo-tools SURVIVOR/1.0.3
#cd $OutPATH
# tried with 100, 10 and 1 caller?(samples)
SURVIVOR merge $OutPATH/sample_files.txt 1000 10 1 1 0 30 $PopulationPath/sample_merge_10.vcf
SURVIVOR merge $OutPATH/sample_files.txt 1000 100 1 1 0 30 $PopulationPath/sample_merge_100.vcf
SURVIVOR merge $OutPATH/sample_files.txt 1000 1 1 1 0 30 $PopulationPath/sample_merge_1.vcf



#TODO: 1.select variants reported by CNVnator, if they have been documented in the Manta population matrix?
# zhiwei94@sens2016007-b9 testSURVIRO_tmp $ SURVIVOR merge
# File with VCF names and paths
# max distance between breakpoints
# Minimum number of supporting caller
# Take the type into account (1==yes, else no)
# Take the strands of SVs into account (1==yes, else no)
# Estimate distance based on the size of SV (1==yes, else no). !!
# Minimum size of SVs to be taken into account.
#cat KA09-155.vcf | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe | wc -l
#10898


####LOGS:
cd /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/SURVIRO
#compare two sigCNVRs reported by two algorithms
bedtools intersect -wao -a /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/SigCNVR_size_changePvalue_withBiomarker.bed -b SigCNVR_size_withBiomarker_Manta.bed | less -S

# test: select PASS (and CNV?)
/sw/bioinfo/GATK/4.1.1.0/bianca/gatk SelectVariants \
-V KA09-155_KA09-0155.vcf \
-O output_KA09-155_KA09-0155.vcf \
----excludeFiltered

bcftools view -f PASS KA09-155_KA09-0155.vcf | less -S
1       869444  MantaDEL:16:0:1:0:0:0
1       1649118 MantaDEL:142:0:0:0:5:0
1_1649118-1_1649170,1_1649597-1_1649649

testStr=AAGCATTTATGTGTACATATGTGTGTTAGCGTGTGCATGCATATGTGTTTGTG
echo ${#testStr}
awk -F '\t' '$1==12 && $2>108971000{print}' /proj/sens2016007/nobackup/VCF-files/Recal_SNP_INDEL_clean_samples_RENAMED.vcf | less -S
# VCF intersect:
head -n3 /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SigCNVR_size_withBiomarker_Manta.bed >  SigCNVR_size_withBiomarker_Manta_head3.bed
cp /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/Association_folder/SigCNVR_size_withBiomarker_Manta.bed ./
#awk -F '\t' '{if($0 ~ /\#/) print; else print $1, $2, $2+1, $0}' sample_merge_100.vcf | less -S
#save header:
awk -F '\t' '{if($0 ~ /\#/) print}' sample_merge_100.vcf > tmp_header_sample_merge_100.vcf
#variant features:
awk '{FS = "\t";OFS = "\t"; if($0 !~ /\#/) print $1,$2,$2+1,$0}' sample_merge_100.vcf > tmp_sample_merge_100.vcf
awk '{FS = "\t";OFS = "\t"; if($0 !~ /\#/) print $1,$2,$2+1,$0}' sample_merge_10.vcf > tmp_sample_merge_10.vcf
awk '{FS = "\t";OFS = "\t"; if($0 !~ /\#/) print $1,$2,$2+1,$0}' sample_merge_1.vcf > tmp_sample_merge_1.vcf
#28Oct new -+5kb
#awk '{FS = "\t";OFS = "\t"; if($0 !~ /\#/) print $1,$2-1000,$2+5000,$0}' sample_merge_1.vcf > tmp_sample_merge_1.vcf
#
module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -b tmp_sample_merge_100.vcf -a SigCNVR_size_withBiomarker_Manta_head3.bed > Intersect_output.vcf
#bedtools intersect -wao -b tmp_sample_merge_10.vcf -a SigCNVR_size_withBiomarker_Manta.bed | less -S
bedtools intersect -wao -b tmp_sample_merge_1.vcf -a SigCNVR_size_withBiomarker_Manta.bed > intersect_merge1_SigManta.bed
#intersect with population matrix: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Population_collapse_bin10.bed
ThePopBed="/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Population_collapse_bin10.bed"
bedtools intersect -wao -b $ThePopBed -a SigCNVR_size_withBiomarker_Manta.bed | awk -v OFS='\t' '{ {Twocount = 0;Onecount=0;Zerocount=0;Threecount=0;Fourcount=0;NAcount=0}; for (i = 9; i <= NF; i++){ if ($i == 4) {Fourcount++};if ($i == 3) {Threecount++};if ($i == ".") {Twocount++}; if ($i == 1) {Onecount++}; if ($i == 0) {Zerocount++};if ($i == "NA") {NAcount++}}; print $1,$2,$3,"4:3:2:1:0:NA", Fourcount ":" Threecount ":" Twocount ":" Onecount ":" Zerocount ":" NAcount }' | less -S
bedtools intersect -wao -b $ThePopBed -a SigCNVR_size_withBiomarker_Manta.bed | awk -v OFS='\t' '{ {Twocount = 0;Onecount=0;Zerocount=0;Threecount=0;Fourcount=0;NAcount=0}; for (i = 9; i <= NF; i++){ if ($i == 4) {Fourcount++};if ($i == 3) {Threecount++};if ($i == ".") {Twocount++}; if ($i == 1) {Onecount++}; if ($i == 0) {Zerocount++};if ($i == "NA") {NAcount++}}; print $1,$2,$3,$4,$6,$7,$8,"4:3:2:1:0:NA", Fourcount ":" Threecount ":" Twocount ":" Onecount ":" Zerocount ":" NAcount }' | bedtools intersect -wao -a stdin -b tmp_sample_merge_1.vcf | less -S
#new
MantBED=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/BinVariants/WinSize10/PasteFolder/Population_collapse_bin10.bed
bedtools intersect -wao -b $MantBED -a /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/SigCNVR_size_changePvalue_withBiomarker.bed > SigCNVnator_inter_MantaBED.bed
awk -v OFS='\t' '{ {Twocount = 0;Onecount=0;Zerocount=0;Threecount=0;Fourcount=0;NAcount=0};for (i = 9; i <= NF; i++){ if ($i == 4) {Fourcount++};if ($i == 3) {Threecount++};if ($i == ".") {Twocount++}; if ($i == 1) {Onecount++}; if ($i == 0) {Zerocount++};if ($i == "NA") {NAcount++}}; print $1,$2,$3,$4,$6,$7,$8,"4:3:2:1:0:NA", Fourcount ":" Threecount ":" Twocount ":" Onecount ":" Zerocount ":" NAcount }'  SigCNVnator_inter_MantaBED.bed > SigCNVnator_inter_MantaBED_Mantafreq.bed
CNVnatorBED=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/PopulationMatrixFolder/changeName_biomarkAssFolder/CNV_1021.collapse.bed
bedtools intersect -wao -b $CNVnatorBED -a SigCNVnator_inter_MantaBED_Mantafreq.bed > SigCNVnator_inter_MantaBED_Mantafreq_CNVnatorBED.bed
awk -v OFS='\t' '{ {Morecount = 0; Twocount = 0; Onecount = 0; Zerocount = 0; NAcount=0}; ;for (i = 13; i <= NF; i++){ if ($i == 2) {Twocount++}; if (0.5 < $i && $i < 2) {Onecount++}; if ($i < 0.5) {Zerocount++}; if ($i > 2) {Morecount++}; if ($i < 0) {NAcount++}}; print $1,$2,$3,$4,$5,$6,$7,$7-$6,$8,$9,$10,$11,$12,$12-$11,"Over2:2:1:0:NA", Morecount ":" Twocount ":" Onecount ":" Zerocount ":" NAcount }'  SigCNVnator_inter_MantaBED_Mantafreq_CNVnatorBED.bed > SigCNVnator_inter_MantaBED_Mantafreq_CNVnatorBED_freq.bed
#compare thwo SigCNVR:
CNVnatorSig=/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/SigCNVR_size_changePvalue_withBiomarker.bed
MantaSig=/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/FormatPopulation/SURVIRO/SigCNVR_size_withBiomarker_Manta.bed
bedtools intersect -a $CNVnatorSig -b $MantaSig
#extend sigCNVRs:
# TheCNVnatorSig="/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/SigCNVR_size_changePvalue_withBiomarker.bed"
# bedtools intersect -wao -b $ThePopBed -a SigCNVR_size_withBiomarker_Manta.bed | awk -v OFS='\t' '{ {Twocount = 0;Onecount=0;Zerocount=0;Threecount=0;Fourcount=0;NAcount=0}; for (i = 9; i <= NF; i++){ if ($i == 4) {Fourcount++};if ($i == 3) {Threecount++};if ($i == ".") {Twocount++}; if ($i == 1) {Onecount++}; if ($i == 0) {Zerocount++};if ($i == "NA") {NAcount++}}; print $1,$2,$3,$4,$6,$7,$8,"4:3:2:1:0:NA", Fourcount ":" Threecount ":" Twocount ":" Onecount ":" Zerocount ":" NAcount }' | bedtools intersect -wao -b stdin -a $TheCNVnatorSig | less -S
# below: only chr11 19 3 6 overlap, the hOSCAR only has 3 individuals
bedtools intersect -wao -b tmp_sample_merge_1.vcf -a /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/SigCNVR_size_changePvalue_withBiomarker.bed > intersect_merge1_SigCNVnator.bed
