# CNVnator_1_FormatBED
## Formattinf CNVnator STDOUT results as BED
## For each individual, CNVnator stdOUT results were saved as .txt file.
## Filter by 0<=q0<=0.5 and converted as BED format
## TXT header line: TXT:$1:CNV_type $2:coordinates(Chr-start-end) $3:CNV_size $4:normalized_RD $5:e-val1 $6:e-val2 $7:e-val3 $8:e-val4 $9:q0
## Converted BED header line: CNVnator_BED:$1:chrom $2:chromStart $3:chromEnd $4:estimatedCN(2*normalized_RD) $5:e-val1(t-test) $6:e-value2() $7:q0, $8:CNV_size
OutPATH=/proj/nobackup/sens2016007/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0
VariantTXT=$(ls -d /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/CNVnator_x*/3_Variant_output/*)

mkdir -p $OutPATH
for CNVnatorCall in ${VariantTXT}
do
  BaseName=${CNVnatorCall##*/}
  awk '{if ($9>=0&&$9<=0.5) print $0}' $CNVnatorCall | awk -v OFS="\t" '{split($2, a, /:|-/); print a[1], a[2]-1, a[3], $4*2, $5, $6, $9, $3}' \
  > ${OutPATH}/${BaseName%.txt}_selected.bed
done
###########################################

# CNVnator_2_FormatPopulationCNVs
## Merge all CNVs in the population as CNVRs
cd $OutPATH
module load bioinfo-tools BEDTools/2.27.1
cat ./*_selected.bed > All.tmp.bed
sort -k1,1 -k2,2n All.tmp.bed > All.tmp.sort.bed
bedtools merge -i All.tmp.sort.bed > All.tmp.sort.merge.bed
sort -k1,1 -k2,2n All.tmp.sort.merge.bed > All.tmp.sort.merge.sort.bed
## Chop the CNVRs by equal size windows
bedtools makewindows -b All.tmp.sort.merge.sort.bed -w 200 > All.pos_200pb_window.bed
ls *selected.bed > Samples.list # list of all bed files
## Map each bed file to fregmented CNVRs, the genotype of a loci is given by the genotype(Copy Number) of the CNVs mapper over 50% of the loci's region.
## For the intersect stdOUT: $7:genotype(copy number) $8:q0-valueX t-test
## If the individual didn't report CNV for a loci, the genotype of the loci for this sample is 2(original diploid genotype)
## If the CNV has low quality (q0X => t-test/ > 0.05/1021), report the genotype as negative value (-1*$7), they were later set as "NA" in GWAS
while read line
do
        echo "bedtools intersect -wao -f 0.5  -a All.pos_200pb_window.bed -b $line | awk '(\$7==\".\"){OFS=\"\\t\";print 2}(\$7!=\".\"&&\$8<=0.05/1021){OFS=\"\\t\";print \$7}(\$7!=\".\"&&\$8>0.05/1021){OFS=\"\\t\";print -1*\$7}' >${line%.bed}_200bp.bed"
done < Samples.list > executeIntersect.sh
# with coordinate: _200bpPOS.bed # add coordinate and some scores:
#_200bpPOS.bed:$1:chrom $2:start $3:end $4:copy_number $5:q0 $6:variant_start $7:variant_end $8:variant_length
###########################################

# CNVnator_3_rmMultipleObservation_onSame_loci
cd $OutPATH
## deDuplicat function: /proj/nobackup/sens2016007/Zhiwei/CNVscripts/dedupForSample.sh
## this step return multiple files for each bed files:
## ${1%.bed}.dedup.bed: is the file containing one genotype for each loci, used for downstream generating CNVs matrix
## ${1%.bed}.DUP.bed: contains the loci that overlap with different CNVs, thus report multiple copy number.
## cat /proj/nobackup/sens2016007/Zhiwei/CNVscripts/dedupForSample.sh
Length=$(wc -l < $1)
DUP=false
Count=0
inLoop=false
Pre_Pos=0 Pre_Start=0 Pre_End=0 Pre_CN=0 Pre_q0=0 Pre_Others=0
# Check if the Others also contain tabs in the string
# name of the output file
echo Length of line in $1: $Length
echo Save dedup file as ${1%.bed}.dedup.bed
echo Duplicate variants: ${1%.bed}.DUP.bed
rm -f ${1%.bed}.DUP.bed
while read Pos Start End CN q0 Others
do
  #echo test $Pos $Start $End $CN
  Count=$((Count + 1))
  if [ "$Count" == $Length ] #&& [ "$DUP" == false ]
  then
    if [ "$DUP" == false ]
    then
      echo -e "$Pre_Pos\t$Pre_Start\t$Pre_End\t$Pre_CN\t$Pre_q0"
      echo -e "$Pos\t$Start\t$End\t$CN\t$q0"
    fi
    if [ "$DUP" == true ]
    then
      echo -e "$Pos\t$Start\t$End\t$CN\t$q0\t$Others" >> ${1%.bed}.tmp.list
      #TODO: only print the first 5 column
      sort -nk5,5 ${1%.bed}.tmp.list | head -n1 | cut -f1-5
      sort -nk5,5 ${1%.bed}.tmp.list >> ${1%.bed}.DUP.bed
    fi
    #break
  elif [ "$Pos" == "$Pre_Pos" ] && [ "$Start" == "$Pre_Start" ] && [ "$End" == "$Pre_End" ]
  then
    if [ "$DUP" == false ]
    then
      echo -e "$Pre_Pos\t$Pre_Start\t$Pre_End\t$Pre_CN\t$Pre_q0\t$Pre_Others" > ${1%.bed}.tmp.list
    fi
    echo -e "$Pos\t$Start\t$End\t$CN\t$q0\t$Others" >> ${1%.bed}.tmp.list
    DUP=true
  elif [ "$DUP" == true ]
    then
      sort -nk5,5 ${1%.bed}.tmp.list | head -n1 | cut -f1-5
      sort -nk5,5 ${1%.bed}.tmp.list >> ${1%.bed}.DUP.bed # change to overlap.bed?
      DUP=false
  elif [ "$Count" != 1 ]
    then
      echo -e "$Pre_Pos\t$Pre_Start\t$Pre_End\t$Pre_CN\t$Pre_q0"
      DUP=false
  fi
  Pre_Pos=$Pos Pre_Start=$Start Pre_End=$End Pre_CN=$CN Pre_q0=$q0 Pre_Others=$Others
done< $1 > ${1%.bed}.dedup.bed
echo number of line after filter: $(wc -l < ${1%.bed}.dedup.bed)
rm -f ${1%.bed}.tmp.list
###########################################

# CNVnator_4_pasteBED as population matrix:
cd $OutPATH
## Only keep the copy number $4 from each sample
for each in *200bpPOS.dedup.bed; do awk '{print $4}' $each > ${each%POS.dedup.bed}.bed; done
cd $OutPATH/subSetFolder/PasteFolder
ls -d $OutPATH/*200bp.bed > Matrix.list
split -l100 Matrix.list --verbose
ls $(pwd)/x* > ListOfGroup.txt

NumLine=$(wc -l < ../../All.pos_200pb_window.bed)
cp $OutPATH/All.pos_200pb_window.bed Population.bed
while read list_1
do
  cat > ${list_1##*/}_Paste.sh <<- EOM
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J ${list_1##*/}_deDUPunfinish
echo "Starting at:"
date
cd $OutPATH/subSetFolder/PasteFolder
NumLine=\$(wc -l < ../../All.pos_200pb_window.bed)
cp ../../All.pos_200pb_window.bed All.pos_200pb_windowMatrix_${list_1##*/}.bed
while read each
  do
    SECONDS=0
    echo Pasting \${each##*/}
    if [ "\$(wc -l < \$each)" != "\$NumLine" ]
    then
      echo "number of lines is different, did you remove overlap variants?"
      break
    fi
    paste All.pos_200pb_windowMatrix_${list_1##*/}.bed \$each > tmpMatrix_${list_1##*/}.bed
    cp tmpMatrix_${list_1##*/}.bed All.pos_200pb_windowMatrix_${list_1##*/}.bed
    duration=\$SECONDS
    echo "Finished, \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."
done < $list_1
rm -f tmpMatrix_${list_1##*/}.bed
EOM
done < /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/subSetFolder/PasteFolder/ListOfGroup.txt
for each in *.sh; do echo "sbatch ${each} > ${each%_*}.id"; done
for each in All.pos_200pb_windowMatrix_xa*.bed; do cut -f 1-3 --complement $each > ${each%.bed}.list; echo "$each done"; done
ls All.pos_200pb_windowMatrix_xa*.bed > ListOfList.txt
# Make a loop to save all the 11 big list to one matrix: Generating pop
cd $OutPATH/subSetFolder/PasteFolder
NumLine=$(wc -l < ../../All.pos_200pb_window.bed)
cp ../../All.pos_200pb_window.bed Population.bed
while read each
do
  SECONDS=0
  echo Pasting ${each##*/}
  if [ "$(wc -l < $each)" != "$NumLine" ]
  then
    echo "number of lines is different, did you remove overlap variants?"
    break
  fi
  paste Population.bed $each > PoptmpMatrix_.bed
  cp PoptmpMatrix_.bed Population.bed
  duration=$SECONDS
  echo "Finished, $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
done < ListOfList.txt
rm PoptmpMatrix_.bed
## Result: /proj/sens2016007/nobackup/Zhiwei/BAManalyses/CNVnatorMain/OptBinSize_q0/subSetFolder/PasteFolder/Population.bed
## Generate header lines, file order in Matrix.list
awk -F" +|/" '{print $1}'  Matrix.list | less -S
awk -OFS='\t' -F" +|/" '{print $NF}'  Matrix.list  > Matrix_basename.list
while read each; do echo ${each%.clean*}; done<Matrix_basename.list > Matrix_shortName.list
## dash line stand for standard input:
printf "#CHROM\nStart\nEnd\n" | cat - Matrix_shortName.list > Matrix_shortName_Pos.list
## Add header to the population.bed file:
paste -sd"\t" Matrix_shortName_Pos.list | cat - Population.bed > Population_withHeader.bed
###########################################

# CNVnator_5_collapseMatrix
## cat:
#!/bin/bash
#SBATCH -A sens2016007
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J collapseMatrix
echo "Starting at:"
date
cd $OutPATH/PopulationMatrixFolder
# define input population matrix formatted by window size:
#inBed="chr19_54399999_54599999.bed"
inBed="Population.bed"
Lines=$(wc -l < $inBed) #! $echo "$Lines" and $echo $Lines give different values
# The population matirx in read into 4 field: Chromosome, StartPosition, EndPosition, Population_genotype(1031 genotype, read as a string)
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
done < $inBed > Population_collapse.bed
## Result:
wc -l Population.bed
##2166502
wc -l Population_collapse.bed
##263831

## Extra remove: /proj/sens2016007/nobackup/work/analyis/ASA_joint_calling_2017_04_03/exlcluded_samples_do_to_genotype_conc.txt
cd $OutPATH/PopulationMatrixFolder
headerLine="$OutPATH/subSetFolder/PasteFolder/Matrix_shortName_Pos.list"
# add header line to the collapse files:
paste -sd"\t" $headerLine | cat - Population_collapse.bed > Population_collapse.rmdup.bed
CollapsedMatrix="Population_collapse.rmdup.bed"
while read each;
do
  echo $each
  # compare to header line:
  #theIndex=$(paste -sd"\t" $headerLine | awk -v each="$each" '{for (i = 4; i <= NF; i++){ if ( $i == each ) print i } }') #not correct! indexing incorrect! add header line to the collapse.bed
  # Compare to population with header.bed:
  theIndex=$(head -n1 $CollapsedMatrix | awk -v each="$each" '{for (i = 4; i <= NF; i++){ if ( $i == each ) print i } }')
  echo "for ${each}, it is in column ${theIndex}, removing..."
  cut -f $theIndex --complement $CollapsedMatrix > ${CollapsedMatrix%.bed}.tmp.bed
  #echo "$theIndex removing the column..." # remove column by awk will keep the columns empty space, need to further format it
  #awk -v OFS='\t' -v theIndex=$theIndex '{$theIndex=""; print $0}' $CollapsedMatrix > ${CollapsedMatrix%.bed}.tmp.bed
  mv ${CollapsedMatrix%.bed}.tmp.bed $CollapsedMatrix
done < /proj/sens2016007/nobackup/work/analyis/ASA_joint_calling_2017_04_03/exlcluded_samples_do_to_genotype_conc.txt
#====> the remove duplicate table: $OutPATH/PopulationMatrixFolder/Population_collapse.rmdup.bed
# Collpase again:
cd $OutPATH/PopulationMatrixFolder
inBed="Population_collapse.rmdup.bed"
Lines=$(wc -l < $inBed) #! $echo "$Lines" and $echo $Lines give different values
# The population matirx in read into 4 field: Chromosome, StartPosition, EndPosition, Population_genotype(1031 genotype, read as a string)
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
done < $inBed > Population_collapse.rmdup.collapsed.bed
wc -l Population_collapse.rmdup.collapsed.bed
#263161
## Before:
wc -l Population_collapse.bed
#263831
###########################################

# CNVnator_6_save_RDS
## Convert sample names:
mkdir changeName_biomarkAssFolder
cp Population_collapse.rmdup.collapsed.bed changeName_biomarkAssFolder/
cd changeName_biomarkAssFolder
# Change the sample names:
# List of the two names Listed in:
nameFile="/proj/sens2016007/nobackup/NSPHS_phenotype_data/WGS_kodnyckel"
while read each id others;
do
  BaseName=${each##*/}
  echo -e "${BaseName%%.*}\t$id"
done < /proj/sens2016007/nobackup/NSPHS_phenotype_data/WGS_kodnyckel > namesDic.txt #TODO: create dictionary?

#head -n 2 namesDic.txt && tail -n +2 namesDic.txt | sort -k1,1 | less -S # sort the table by the first column, same order as the given excel file
head -n1 Population_collapse.rmdup.collapsed.bed | less -S
# Get the header line
head -n1 Population_collapse.rmdup.collapsed.bed > rmdup_header.txt
# Convert to other name:
for each in $(cat rmdup_header.txt);
do
  #echo "Filename $each"
  while read Filename Id;
  do
    if [ "$each" = "$Filename" ]
    then
      echo "$Id"
    fi
  done < namesDic.txt
done > rmdup_header_changeName.txt
printf "CHROM\nStart\nEnd\n" | cat - rmdup_header_changeName.txt > rmdup_header_changeName_Pos.txt
## Add header to the population.bed file:
paste -sd"\t" rmdup_header_changeName_Pos.txt | cat - Population_collapse.rmdup.collapsed.bed > Population_collapse.rmdup.collapsed.rename.bed   # contain both name
## Remove the second row and save it as new file:
sed '2d' Population_collapse.rmdup.collapsed.rename.bed > CNV_1021.collapse.bed
## In R, save CNV_1021.collapse.bed as rds file:
CNVbed <- read.table("CNV_1021.collapse.bed", header = TRUE, check.names=FALSE, sep = "\t") # the "-" is read as ".", use check.names
## Save as RDS:
saveRDS(CNVbed,"CNVbed.rds")
