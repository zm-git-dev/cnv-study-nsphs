#!Rscript
library(data.table)
library(stringr)
library(magrittr)

is_not_empty <- function(dt){
  nrow(dt) > 0
}

filenames <- "!{bed_files}" %>%
  str_split(" ") %>%
  unlist()

sample_ids <- str_split_fixed(filenames, "\\.", n=2)[,1]
names(filenames) <- sample_ids

read_bed <- function(filename) {
  bed <- fread(filename)
  if (nrow(bed) > 0) as.numeric(bed[,4]) else NULL
}

genotypes <- filenames %>%
  lapply(read_bed) %>%
  as.data.table()

first_table <- fread(filenames[1])
genotypes <- first_table[, .(chr = V1, start = V2, end=V3)] %>%
  cbind(genotypes)

write.table(genotypes, "cnv_matrix.txt", sep="\t", quote=F, row.names=F)