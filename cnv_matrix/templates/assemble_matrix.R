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

genotypes_long <- filenames %>%
  lapply(fread) %>%
  Filter(is_not_empty, .) %>%
  lapply(setnames, 1:4, c("chr", "start", "end", "cn")) %>%
  rbindlist(idcol="sample")

genotypes <- dcast(chr + start + end ~ sample, value.vars="cn")

write.table(genotypes, "cnv_matrix.txt", sep="\t", quote=F, row.names=F)