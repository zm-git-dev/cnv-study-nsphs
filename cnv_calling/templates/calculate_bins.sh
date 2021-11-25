#!/bin/bash
# load the ROOT library
module load ROOT/6.06.08 bioinfo-tools CNVnator

bin_sizes=(70 85 100 150 200 250)
chrom=$(seq 1 22)

for bin_size in $bin_sizes
do
  cnvnator -root "!{root}" -his "$bin_size" -chrom $chrom -d "!{reference}"
  # CNVnator outputs a log message from which we have to extract the mean RD and its SD.
  echo "$(cnvnator -root "!{root}" -stat "$bin_size" -chrom $chrom | perl -ne 'print "$1 $2\n" if /^Average RD per bin \(1-22\) is (\d+\.\d+) \+- (\d+\.\d+) \(after GC correction\)/') $bin_size" >> stats.txt
done

max_ratio=5
best_size=-1
while read -r mean sd size; do
  ratio=$(echo "$mean / $sd" | bc -l)
  echo $ratio $size
  if [[ $ratio -le $max_ratio && $ratio -ge 4 ]]; then
    max_ratio=$ratio
    best_size=$size
  fi
done < stats.txt
