#!/bin/bash
module load ROOT/6.06.08 bioinfo-tools CNVnator
chrom=$(seq 1 22)
bin_size=$(echo "!{bin_size}" | xargs)
cnvnator -root !{root} -partition "${bin_size}" -chrom $chrom