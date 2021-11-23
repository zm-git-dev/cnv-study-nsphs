The two individuals with pacBio data is
KA06_0017 and KA06_0031
# the CNVnator and Manta variants calling is documented in 1_VariantsCalling.sh in the FinalReportCode
# List of the two names Listed in:
nameFile="/proj/sens2016007/nobackup/NSPHS_phenotype_data/WGS_kodnyckel"
/proj/sens2016007/delivery/P2656/P2656_342/03-BAM/P2656_342.clean.dedup.bam     KA06-0017       P2656_342_P2656_342
/proj/sens2016007/delivery/P2656/P2656_312/03-BAM/P2656_312.clean.dedup.bam     KA06-0031       P2656_312_P2656_312
# which of the subgroup this two individuals are in ?

# get CNVnator info from the population CNV matrix or individual results?
#CNVnator:
find /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0 -name "P2656_342*"
find /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0 -name "P2656_312*"

# CNVnator raw output:
find /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/3_Variant_output/* -name "P2656_342*"
/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xaa/3_Variant_output/P2656_342.clean.dedup.variant.txt
find /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/3_Variant_output/* -name "P2656_342*"
/proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xaa/3_Variant_output/P2656_342.clean.dedup.variant.txt
less -S /sw/apps/bioinfo/CNVnator/0.3.3/bianca/cnvnator2VCF.pl
ReferenceGenome=/proj/sens2016007/nobackup/Reference/human_g1k_v37.fasta
# directory with individual reference fasta files:
/proj/sens2016007/nobackup/Reference
# test:
module load bioinfo-tools SURVIVOR
# tried to format as VCF by survivor AND CNVnator scripts
cp /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_xaa/3_Variant_output/P2656_342.clean.dedup.variant.txt ./
/sw/apps/bioinfo/CNVnator/0.3.3/bianca/cnvnator2VCF.pl -prefix study1 P2656_342.clean.dedup.variant.txt \
/proj/sens2016007/nobackup/Reference > testCNVnator_bedtovcf.vcf

# e.g.:
P2656_312.clean.dedup.variant_selected.bed # filter by 0 <= q0 <= 0.5
#$1:chrom $2:chromStart $3:chromEnd $4:estimatedCN(2*normalized_RD) $5:e-val1(t-test) $6:e-value2() $7:q0, $8:CNV_size
#P2656_312.clean.dedup.variant_selected_200bp.bed # chopped into 200 bp window
#P2656_312.clean.dedup.variant_selected_200bpPOS.bed # add coordinate and some scores:
#$1:chrom $2:start $3:end $4:copy_number $5:q0 $6:variant_start $7:variant_end $8:variant_length
#P2656_312.clean.dedup.variant_selected_200bpPOS.dedup.bed # remove windows with multiple record due to overlap
#"$Pre_Pos\t$Pre_Start\t$Pre_End\t$Pre_CN\t$Pre_q0"
#P2656_312.clean.dedup.variant_selected_200bpPOS.DUP.bed # the overlapping windows, not all samples have this file.


cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/TwoPacBioSamples
# SelectedBED: select by p0<0.5 TODO: t-test select and format results of two softwares at individual level:
#selectedBED: #$1:chrom $2:chromStart $3:chromEnd $4:estimatedCN(2*normalized_RD) $5:e-val1(t-test) $6:e-value2() $7:q0, $8:CNV_size
zhiwei94@sens2016007-b9 TwoPacBioSamples $ cp /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/P2656_342.clean.dedup.variant_selected.bed ./
zhiwei94@sens2016007-b9 TwoPacBioSamples $ cp /proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/P2656_312.clean.dedup.variant_selected.bed ./
CNVnatorSample1="P2656_342.clean.dedup.variant_selected.bed"
CNVnatorSample2="P2656_312.clean.dedup.variant_selected.bed"
#TODO: CNVnator lable low t-test
# !!label (t-test p-val/ > 0.05/1021) features as negaive CN
module load bioinfo-tools BEDTools/2.27.1
# Old method of formatting population observation: ==> use both CNVnator and Manta to have the 2 individuals 2 algorithms table?
#cat ./*_selected.bed > All.tmp.bed
#sort -k1,1 -k2,2n All.tmp.bed > All.tmp.sort.bed
#bedtools merge -i All.tmp.sort.bed > All.tmp.sort.merge.bed
#sort -k1,1 -k2,2n All.tmp.sort.merge.bed > All.tmp.sort.merge.sort.bed
#bedtools intersect -wao -f 0.5  -a All.tmp.sort.merge.sort.bed -b $Sample1 | awk '($7=="."){OFS="\t";print $1, $2, $3, $3-$2, $4, $5, $6, $6-$5, 2}($7!="."&&$8<=0.05/1021){OFS="\t";print $1, $2, $3, $3-$2, $4, $5, $6, $6-$5, $7}($7!="."&&$8>0.05/1021){OFS="\t";print $1, $2, $3, $3-$2, $4, $5, $6, $6-$5, -1*$7}' \
#> ${Sample1%.bed}_Formatted.bed
# NEW: t-test filter at individual level: two samples: adjusted $5:p-value (0.05/2)????!
for Call in $CNVnatorSample1 $CNVnatorSample2
do
  FileName=CNVnator_${Call%%.*}.PASS_T_test.bed
  awk '($5 <=0.05/2){OFS="\t"; print $1,$2,$3,$3-$2,$4,"PASS"}($5 >0.05/2){OFS="\t"; print $1,$2,$3,$3-$2,$4,"T_test"}' $Call > $FileName
done
awk '($5 <=0.05/2){OFS="\t"; print $1,$2,$3,$3-$2,$4,"PASS"}($5 >0.05/2){OFS="\t"; print $1,$2,$3,$3-$2,$4,"T_test"}' $CNVnatorSample1 > CNVnator_${BaseName%%.}.PASS_Ploidy.bed

#Manta1.5: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main
# For Manta results, keep PASS and !Ploidy! variants.
# header diploidSV.vcf.gz: #CHROM_A        START_A END_A   CHROM_B START_B END_B   ID      QUAL    STRAND_A        STRAND_B        TYPE    FILTER  INFO    FORMAT  P2656_34
#find /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main -name "P2656_342*"
MantaSample1="/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaa/IntermediaOutput/P2656_342_Folder/results/variants/diploidSV.vcf.gz"
#find /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main -name "P2656_312*"
MantaSample2="/proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaa/IntermediaOutput/P2656_312_Folder/results/variants/diploidSV.vcf.gz"

#cp /proj/sens2016007/nobackup/Zhiwei/BAManalyses/Manta1_5Main/Manta_xaa/IntermediaOutput/P2656_342_Folder/results/variants/diploidSV.vcf.gz ./
zcat diploidSV.vcf.gz | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe |  awk '(($11=="DUP"||$11=="DEL")&&($12=="Ploidy"||$12=="PASS")&&$1==$4&&$6>$2){OFS="\t";print $1,$2,$6,$6-$2,$11,$12,$14,$15}(($11=="DUP"||$11=="DEL")&&$12=="PASS"&&$1==$4&&$6<$2){OFS="\t";print $1,$5,$3,$3-$5,$7,$11,$12,$14,$15}' | less -S

for Call in $MantaSample1 $MantaSample2
do
  BaseName=$(basename $(echo ${Call} | sed 's/.*\(IntermediaOutput\)/\1/g' | cut -d'/' -f-2) )
  echo Saving to Manta_${BaseName%_Folder}.PASS_Ploidy.bed
  # header: Chr, Start, End, Size, SVtype, Filter, FORMAT
  zcat $Call | /proj/sens2016007/nobackup/Tools/svtools/vcfToBedpe | awk '(($11=="DUP"||$11=="DEL")&&($12=="Ploidy"||$12=="PASS")&&$1==$4&&$6>$2){OFS="\t";print $1,$2,$6,$6-$2,$11,$12,$14,$15}(($11=="DUP"||$11=="DEL")&&($12=="Ploidy"||$12=="PASS")&&$1==$4&&$6<$2){OFS="\t";print $1,$5,$3,$3-$5,$11,$12,$14,$15}' > Manta_${BaseName%_Folder}.PASS_Ploidy.bed
done

cd KA06_0017_P2656_342_Folder
# Manta filter by size?
# study the distribution of variant size of Manta: the top 5% lines:
sort -nrk4 Manta_P2656_342.PASS_Ploidy.bed | head -n235 | less -S
# filter Manta Size$4 by 2Mb: !!
awk '{if ($4 <= 2000000) {print}}' Manta_P2656_342.PASS_Ploidy.bed > Manta_P2656_342.PASS_Ploidy_Size2M.bed
#cat *.bed > All.tmp_bed
#sort -k1,1 -k2,2n All.tmp_bed > All.tmp.sort.bed
module load bioinfo-tools BEDTools/2.27.1
cat CNVnator_P2656_342.PASS_T_test.bed Manta_P2656_342.PASS_Ploidy_Size2M.bed | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3 }' | bedtools merge | awk '{OFS="\t"; print $1,$2,$3,$3-$2}' > CNVR.bed
#awk '{OFS="\t"; print $1, $2, $3}' All.tmp.sort.bed | bedtools merge | awk '{OFS="\t"; print $1,$2,$3,$3-$2}' > All.tmp.sort.merge.bed
bedtools intersect  -wao -a CNVR.bed -b CNVnator_P2656_342.PASS_T_test.bed | bedtools intersect -wao -a stdin -b Manta_P2656_342.PASS_Ploidy_Size2M.bed > CNVR_CNVnator_Manta_P2656_342.bed
# get SignificantCNVR:
SigCNVR="/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed"
bedtools intersect -wao -a $SigCNVR -b CNVR_CNVnator_Manta_P2656_342.bed | less -S # then saved as SigCNVR...

mkdir KA06_0031_P2656_312_Folder
cp CNVnator_P2656_312.PASS_T_test.bed ./KA06_0031_P2656_312_Folder/
cp Manta_P2656_312.PASS_Ploidy.bed ./KA06_0031_P2656_312_Folder/
cd KA06_0031_P2656_312_Folder
CNVnatorFile="CNVnator_P2656_312.PASS_T_test.bed"
MantaFile="Manta_P2656_312.PASS_Ploidy.bed"
# Filter Manta by size:
#sort -nrk4 Manta_P2656_312.PASS_Ploidy.bed | head -n220 | less -S #220 estimated by 5% wc -l
awk '{if ($4 <= 2000000) {print}}' $MantaFile > ${MantaFile%.bed}_Size2M.bed
cat $CNVnatorFile ${MantaFile%.bed}_Size2M.bed | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3 }' | bedtools merge | awk '{OFS="\t"; print $1,$2,$3,$3-$2}' > CNVR.bed
bedtools intersect  -wao -a CNVR.bed -b $CNVnatorFile | bedtools intersect -wao -a stdin -b ${MantaFile%.bed}_Size2M.bed > CNVR_CNVnator_${MantaFile%%.*}.bed
# get SignificantCNVR:
SigCNVR="/proj/nobackup/sens2016007/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/SignificantCNPs_merged_size.bed"
bedtools intersect -wao -a $SigCNVR -b CNVR_CNVnator_${MantaFile%%.*}.bed | less -S
bedtools intersect -wao -a $SigCNVR -b CNVR_CNVnator_${MantaFile%%.*}.bed > SigCNVR_localCNVR_both_${MantaFile%%.*}.bed

cd TwoSampleIntersect/
cat ../KA06_0017_P2656_342_Folder/CNVnator_P2656_342.PASS_T_test.bed \
../KA06_0017_P2656_342_Folder/Manta_P2656_342.PASS_Ploidy_Size2M.bed \
../KA06_0031_P2656_312_Folder/CNVnator_P2656_312.PASS_T_test.bed \
../KA06_0031_P2656_312_Folder/Manta_P2656_312.PASS_Ploidy_Size2M.bed | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3 }' | bedtools merge | awk '{OFS="\t"; print $1,$2,$3,$3-$2}' > Both_CNVR.bed

bedtools intersect -wao -a Both_CNVR.bed -b ../KA06_0017_P2656_342_Folder/CNVnator_P2656_342.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0017_P2656_342_Folder/Manta_P2656_342.PASS_Ploidy_Size2M.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/CNVnator_P2656_312.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/Manta_P2656_312.PASS_Ploidy_Size2M.bed | less -S
# detect by both algorithms and both individuals:
# targeted columns: $5 $12 21 28
bedtools intersect -wao -a Both_CNVR.bed -b ../KA06_0017_P2656_342_Folder/CNVnator_P2656_342.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0017_P2656_342_Folder/Manta_P2656_342.PASS_Ploidy_Size2M.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/CNVnator_P2656_312.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/Manta_P2656_312.PASS_Ploidy_Size2M.bed | \
awk '{if (($5 != ".")&&($12 != ".")&&($21 != ".")&&($28 != ".")) {print}}' | less -S
# more stringent Manta size filter?
# filter Plodity to have less nested vatiants?
bedtools intersect -wao -a Both_CNVR.bed -b ../KA06_0017_P2656_342_Folder/CNVnator_P2656_342.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0017_P2656_342_Folder/Manta_P2656_342.PASS_Ploidy_Size2M.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/CNVnator_P2656_312.PASS_T_test.bed | \
bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/Manta_P2656_312.PASS_Ploidy_Size2M.bed | \
awk '{if (($5 != ".")&&($12 != ".")&&($21 != ".")&&($28 != ".")&&($17 != "Ploidy")&&($33 != "Ploidy")) {print}}' | less -S
# the size still vary a lot!
###TODO: overlapping regions set as NA and check the percentage of those regions in all variable regions.
less -S CNVR_CNVnator_Manta_P2656_342.bed
621431 # current filter threshold 2Mb, change? just use GQ to filter?: 2 600k+ and one 2k variants reported by Manta
209935006 #  # in document <15 as low GQ, report DUP(742bp) and DEL(1069bp) in a 1070bp region: both deletion and duplication around the same regions, with high quality!
63698905 # on chr11: same, DUP overlap with DEL, around same size. 
93019920 # CNVnator suggests het DEL, Manta: higher GQ for het DEL in the later two results
8558480 # on chr12, CNVnator suggests homo DEL, Manta: the later 3 results overlap, higher for homo DEL
awk '{print $1, $2, $3, $4}' CNVR_CNVnator_Manta_P2656_342.bed | sort -rnk4,4 | uniq | less -S

# check the two individuals have the sigCNVRs:
changeP_SigCNVR="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
changeP_SigCNVR_biomarker="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue_withBiomarker.bed"
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/TwoPacBioSamples/KA06_0017_P2656_342_Folder
KA06_0017_P2656_342_Folder $ bedtools intersect -wao -a $changeP_SigCNVR_biomarker -b CNVnator_P2656_342.PASS_T_test.bed | bedtools intersect -wao -a stdin -b Manta_P2656_342.PASS_Ploidy_Size2M.bed | bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/CNVnator_P2656_312.PASS_T_test.bed | bedtools intersect -wao -a stdin -b ../KA06_0031_P2656_312_Folder/Manta_P2656_312.PASS_Ploidy_Size2M.bed > ../TwoSampleIntersect/Both_changeP_interestct.bed
less -S ../TwoSampleIntersect/Both_changeP_interestct.bed
cd ../TwoPacBioSamples/
#awk '{OFS="\t"; print $1,$2,$3,$4,$9,$10,$11,$16,$17,$18,$25,$26,$27,$32,$33,$34,$5}' Both_changeP_interestct.bed > Compact_Both_changeP_interestct.bed
# format the genotype(GT) info reported by Manta:
cat Both_changeP_interestct.bed | awk '{$20 =substr($20, 1, 3);$36 =substr($36, 1, 3); print $20,$36}' | less -S
# header: CNVR(Chr,Star,End,Size); KA06-0017(CNVnator:Size,CN,FT; Manta:Size,CN,GT,FT;); KA06-0031(CNVnator:Size,CN,FT; Manta:Size,CN,GT,FT;); Biomarker:numberOfWin
awk '{OFS="\t";$20 =substr($20, 1, 3);$36 =substr($36, 1, 3); print $1,$2,$3,$4,$9,$10,$11,$16,$17,$20,$18,$25,$26,$27,$32,$33,$36,$34,$5}' Both_changeP_interestct.bed > Compact_Both_changeP_interestct.bed
cp Compact_Both_changeP_interestct.bed /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/SummerProject/SigCNVR_changeP_folder/.
## make deletions and duplications table for this two individuals? or just seperate DEL and DUP for Manta population results.
