module load ROOT/6.06.08 bioinfo-tools CNVnator
chrom=$(seq 1 22)
echo Reading BAM file $bam
cnvnator -root out -chrom $chrom -tree $bam -unique