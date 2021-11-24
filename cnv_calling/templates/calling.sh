module load ROOT/6.06.08 bioinfo-tools CNVnator
variants_file="$(basename !{root} .root)_variants.txt"
cnvnator -root !{root} -call !{bin_size} > $variants_file