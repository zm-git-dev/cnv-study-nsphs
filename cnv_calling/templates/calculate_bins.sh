#!/bin/bash
# load the ROOT library
module load ROOT/6.06.08 bioinfo-tools CNVnator

bin_sizes=(70 85 100 150 200 250)
chrom=$(seq 1 22)

max_ratio=5
best_size='-1'
for bin_size in bin_sizes
do
  cnvnator -root "!{root}" -his "$bin_size" -chrom $chrom -d "!{referece}"
  # CNVnator outputs a log message from which we have to extract the mean RD and its SD.
  bin_stats=($(cnvnator -root "!{root}" -stat "$bin_size" -chrom $chrom | perl -ne 'print "$1\t$2\n" if /^Average RD per bin \(1-22\) is (\d+\.\d+) \+- (\d+\.\d+) \(after GC correction\)/'))
  let "ratio=${bin_stats[1]}/${bin_stats[2]}"
  # The optimum ratio is between 4 and 5. The smaller, the better.
  # Therefore, the smallest found ratio between 4 and 5 is the upper bound for the optimum
  if [[ "$ratio" -le "$max_ratio" && "$ratio" -ge 4 ]]; then
    max_ratio=$ratio
    best_size=$bin_size
  fi
done
echo $best_size