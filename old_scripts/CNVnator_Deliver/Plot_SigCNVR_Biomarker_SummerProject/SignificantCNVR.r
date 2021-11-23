#22 JUL
#mkdir /proj/sens2016007/nobackup/Zhiwei/SummerProject
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
# in Bash
less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/Sig_CNVR_info/Intrsect.bed

# in R

#rm(list=ls())
# read pea_3 data:
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
# ls()
AgeSex <- read.csv(file="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv", header=TRUE, sep=",")
Intersect_bed <- read.table("/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/GWAS/Sig_CNVR_info/Intrsect.bed", header = FALSE, check.names=FALSE, sep = "\t")

# !remove the last 2 columns from bedtools intersect overlap region SIZE, add header to the data frame.
Intersect_bed <- subset(Intersect_bed, select = -c(889,890))
# show the first 18 columns
Intersect_bed[1:2,1:18]
Intersect_bed <- subset(Intersect_bed, select = -c(7,8,9,14,15,16))

# Read header ==> read the first line
#less -S /proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt
CNVbed <- read.table("/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt", sep="\t", col.names=TRUE, row.names=FALSE)
#Error in read.table("/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt",  :
#  more columns than column names ?????
#return a string: readLines("/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt", n = 1)
# default comment.char = "#"
# only used for getting the sample names!
CNVbed <- read.table("/proj/sens2016007/nobackup/wharf/zhiwei94/zhiwei94-sens2016007/tables/tables/CNVnatorTables/CNdiatribute/CNP_filtr.txt", sep="\t", header = TRUE, comment.char = "", check.names = FALSE)
SampleName <- colnames(CNVbed)
rm(CNVbed)
# column names
c(SampleName[1:3], "Size", "GLMresult", "Biomarker", "CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[4:length(SampleName)])

# add header line to Intersect_bed
colnames(Intersect_bed) <- c(SampleName[1:3], "Size", "GLMresult", "Biomarker", "CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[4:length(SampleName)])
unique(Intersect_bed[,6:10]) # inclusing biomarker
unique(Intersect_bed[,7:10]) # only count CNVRs
#23 Jul
# Group this table, for copy number histogram and CN~Biomarker distribution, excluding biomarker information
deduped.Intersect_bed <- unique( Intersect_bed[ , c(1:4, 7:ncol(Intersect_bed)) ] )
# Select the Significant CNVR:
Sig_CNVR <- unique( deduped.Intersect_bed[ , c(5:8) ] )

# select by multiple columns(CNVR) be the same as loop in Sig_CNVR
##! error message:
# Within the loop:
i <- 1
loop_CNVR <- Sig_CNVR[i,]
loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
# test selecting the window in the same CNVR
# ggplot histogram for the copy number genotype from the same CNVR, from $9, only the copy number info
loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]
CopyNumber <- c(t(loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]))
CN.df <- data.frame(CopyNumber)
ggplot(data= CN.df , aes(x= CopyNumber)) +
   geom_histogram()
ggsave(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, ".png"))
library("ggplot2")
> ggplot(data= CN.df , aes(x= CopyNumber)) +
+   geom_histogram()
ggsave("test.png")
# CN=2 has the highest frequency, but hard to visualize the other CNs
# 24 JUL, define function:
PlotCopyNumber <- function(loop_CNVR, deduped.Intersect_bed) {
  # select window in the same CNVR
  loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
  # ggplot histogram for the copy number genotype from the same CNVR, from $9, only the copy number info
  loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]
  CopyNumber <- c(t(loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]))
  CN.df <- data.frame(CopyNumber)
  ggplot(data= CN.df , aes(x= CopyNumber)) +
     geom_histogram()
  ggsave(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, ".png"))
}


# Copy Number distribution plot, !only CN distribution!
i <- 1
loop_CNVR <- Sig_CNVR[i,]
PlotCopyNumber(loop_CNVR, deduped.Intersect_bed)
# Plot CN distribution for all significant 30 CNVRs
for ( i in 1:nrow(Sig_CNVR) ) {
  loop_CNVR <- Sig_CNVR[i,]
  PlotCopyNumber(loop_CNVR, deduped.Intersect_bed)
}

# 24 JUL 26 JUL, TODO: change the index to the variable to columns names, maybe easiler to read?
# Individual with both DNA and biomarker data:
BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
# Sort by sample names' order
BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
# test select the associated biomarker in the same CNVR
i <- 1 # index for CNVR
# select one CNVR
loop_CNVR <- Sig_CNVR[i,]
# select window in the same CNVR
loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
# select associate biomarker in the same CNVR
loop_Biomarker <- unique(Intersect_bed[ Intersect_bed[,7] == as.numeric(loop_CNVR[1]) & Intersect_bed[,8] == as.numeric(loop_CNVR[2]) & Intersect_bed[,9] == as.numeric(loop_CNVR[3]), ][,6])
# !TODO: make loop to go through all biomarker in loop_Biomarker
# TODO: keep the original regression result in the intersect_bed file?
# generate data frame for making linear regression plot
# foo loop: Index <- 1
IndexBiomarker <- 1 # index for biomarker: length(loop_Biomarker)
# double in the case of loop_Biomarker contains more than 1 value, can be index correctly
One_Biomarker <- loop_Biomarker[1]
# Select biomarker info:
# !no longer a dataframe if without data.frame()
One_Biomarker_value <- data.frame(BioMarker.df[, One_Biomarker])
# add biomarker name to column name
colnames(One_Biomarker_value) <- One_Biomarker

# merge CN and Biomarker_value, for linear regression analysis:
# TODO: transpose the One_Biomarker_value table # see next section
# 7 Aug: plot biomarker with CN, with plot or boxplot function
# select the CNV table WITHOUT coordinate info of window and CNVRs
loop_deduped.Intersect.noCoor <- loop_deduped.Intersect[,-(1:8)]
loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect.noCoor[,order(colnames(loop_deduped.Intersect.noCoor))])
# thr sorted and transposed table should have same same rows as dim(One_Biomarker_value)
#png('test_CNV_biomarker.png')
# By saving into pdf, it will automatically append to the sme file
####### R, plot CN vs Biomarker, dot plot
pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".Test.pdf"))
for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
  plot(cbind(data.frame(loop_deduped.Intersect.sort.t[,value]), list(One_Biomarker_value)), type = "p",
  main = (paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "_Biomarker:", as.character(One_Biomarker) )),
  xlab = (paste0( loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_WindowSize:", loop_deduped.Intersect[value,]$Size)))
}
dev.off()
####### ggplot, boxplot
library(ggplot2)
value <- 1
# Use aes_string instead of aes for indexing the columns for X and Y
one.CN_Biomarker.df <- data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
one.CN_Biomarker.df$CN <- as.factor(one.CN_Biomarker.df$CN)
# ! below: have to use the column direct name to index the biomarker values!
#p <- ggplot(one.CN_Biomarker.df, aes(x=CN, y=ONC2_194_FCRLB)) +
p <- ggplot(one.CN_Biomarker.df, aes_string(x=names(one.CN_Biomarker.df)[1], y=names(one.CN_Biomarker.df)[2])) +
  geom_boxplot() +
  ggtitle(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "_Biomarker:", as.character(One_Biomarker) )) +
  labs(x = paste0( loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_WindowSize:", loop_deduped.Intersect[value,]$Size, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t) ))
ggsave("myGGplotBoxTest2.png", plot = p)
####### move into loop to format on PDF for one CNVR:
pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".rounded.boxplot.pdf"))
for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
  one.CN_Biomarker.df <- data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
  one.CN_Biomarker.df$CN <- as.factor(one.CN_Biomarker.df$CN)
  p <- ggplot(one.CN_Biomarker.df, aes_string(x=names(one.CN_Biomarker.df)[1], y=names(one.CN_Biomarker.df)[2])) +
    geom_boxplot() +
    ggtitle(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "_Biomarker:", as.character(One_Biomarker) )) +
    labs(x = paste0( loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_WindowSize:", loop_deduped.Intersect[value,]$Size, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t) ))
  #ggsave("myGGplotBoxTest2.png", plot = p)
  print(p)
}
dev.off()


# make function to process each Sig_CNVR input function to make dot plot
PlotCN_Biomarker <- function(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3) {
  BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
  # Sort by sample names' order
  BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
  # select window in the same CNVR
  loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
  print(paste0("Number of CNVR-window: ",(dim(loop_deduped.Intersect)[1]) ))
  #print(dim(loop_deduped.Intersect))
  # select associate biomarker in the same CNVR
  loop_Biomarker <- unique(Intersect_bed[ Intersect_bed[,7] == as.numeric(loop_CNVR[1]) & Intersect_bed[,8] == as.numeric(loop_CNVR[2]) & Intersect_bed[,9] == as.numeric(loop_CNVR[3]), ][,6])
  for (BiomarkerIndex in 1:length(loop_Biomarker)) {
    One_Biomarker <- loop_Biomarker[BiomarkerIndex]
    One_Biomarker_value <- data.frame(BioMarker.df[, One_Biomarker])
    colnames(One_Biomarker_value) <- One_Biomarker
    #loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect[,-(1:8)][order(colnames(loop_deduped.Intersect[,-(1:8)]))])
    loop_deduped.Intersect.noCoor <- loop_deduped.Intersect[,-(1:8)]
    loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect.noCoor[,order(colnames(loop_deduped.Intersect.noCoor))])
    # round copy number: if not rounded, will be scattered for not CN=2 regions
    loop_deduped.Intersect.sort.t <- round(loop_deduped.Intersect.sort.t)
    #pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".pdf"))
    pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".rounded.pdf"))
    for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
      plot(cbind(data.frame(loop_deduped.Intersect.sort.t[,value]), list(One_Biomarker_value)), type = "p",
      main = (paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "_Biomarker:", as.character(One_Biomarker) )),
      xlab = (paste0( loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_WindowSize:", loop_deduped.Intersect[value,]$Size, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t)  )))
    }
    dev.off()
  }
}

i <- 2
loop_CNVR <- Sig_CNVR[i,]
PlotCN_Biomarker(loop_CNVR, deduped.Intersect_bed, Intersect_bed)

for ( i in 1:nrow(Sig_CNVR) ) {
  loop_CNVR <- Sig_CNVR[i,]
  PlotCN_Biomarker(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3)
}

###Report script:
#######Breakpoint visualiztion:
#Target: loop_deduped.Intersect
library(ggplot2)

Cor.CV.table = data.frame()
for (rowIndex in 1:nrow(loop_deduped.Intersect)) {
  Cor.CV.df <- data.frame( rep(loop_deduped.Intersect[rowIndex,2], ncol(loop_deduped.Intersect)-8 ), t(loop_deduped.Intersect[rowIndex, 9:ncol(loop_deduped.Intersect)]) )
  colnames(Cor.CV.df) <- c("Coor", "CN")
  rownames(Cor.CV.df) <-  paste(rownames(Cor.CV.df), "Win", rowIndex, sep = "_")
  Cor.CV.table <- rbind(Cor.CV.table, Cor.CV.df)
}
Cor.CV.table$Coor <- as.factor(Cor.CV.table$Coor)
p <- ggplot(Cor.CV.table, aes_string(x=names(Cor.CV.table)[1], y=names(Cor.CV.table)[2])) +
  geom_boxplot() +
  ggtitle(paste0("Copy number boxplot for CNVR:", loop_deduped.Intersect$CNVRchr[1], ":", loop_deduped.Intersect$CNVRstart[1], "-", loop_deduped.Intersect$CNVRend[1], "_Size:", loop_deduped.Intersect$CNVRsize[1], "bp," "_noWin:",  nrow(loop_deduped.Intersect))) +
  labs(x = "Coordinate", y = "Copy Number" )
ggsave("testBreakPointCN.png")

PlotCopyNumber_boxplot <- function (loop_deduped.Intersect) {
  Cor.CV.table = data.frame()
  for (rowIndex in 1:nrow(loop_deduped.Intersect)) {
    # only pick the starting position of a window
    Cor.CV.df <- data.frame( rep(loop_deduped.Intersect[rowIndex,2], ncol(loop_deduped.Intersect)-8 ), t(loop_deduped.Intersect[rowIndex, 9:ncol(loop_deduped.Intersect)]) )
    colnames(Cor.CV.df) <- c("Coor", "CN")
    rownames(Cor.CV.df) <-  paste(rownames(Cor.CV.df), "Win", rowIndex, sep = "_")
    Cor.CV.table <- rbind(Cor.CV.table, Cor.CV.df)
  }
  Cor.CV.table$Coor <- as.factor(Cor.CV.table$Coor)
  p <- ggplot(Cor.CV.table, aes_string(x=names(Cor.CV.table)[1], y=names(Cor.CV.table)[2])) +
    geom_boxplot() +
    ggtitle(paste0("CN boxplot for CNVR:", loop_deduped.Intersect$CNVRchr[1], ":", loop_deduped.Intersect$CNVRstart[1], "-", loop_deduped.Intersect$CNVRend[1], "_Size:", loop_deduped.Intersect$CNVRsize[1], "bp", "_NoWin:",  nrow(loop_deduped.Intersect))) +
    labs(x = "Coordinate", y = "Copy Number", subtitle = "* the width of each box is not scaled to the window's size" )
  print(p)
}

#######Copy number histogram subfunction:
PlotCopyNumberSub <- function(loop_CNVR, deduped.Intersect_bed ) {
  # select window in the same CNVR
  loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
  # ggplot histogram for the copy number genotype from the same CNVR, from $9, only the copy number info
  #loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]
  CopyNumber <- c(t(loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]))
  CN.df <- data.frame(CopyNumber)
  CN.p <- ggplot(data= CN.df , aes(x= CopyNumber)) +
     geom_histogram() +
     ggtitle(paste0("CN histogram for CNVR:", loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "bp", "_NoWin:", nrow(loop_deduped.Intersect) )) +
     labs(x = "Copy Number")
  print(CN.p)
}

########add text to the CN_Biomarker histogram:

#give.n <- function(x){
#  return(c(y = median(x)*1.05, label = length(x)))
#  # experiment with the multiplier to find the perfect position
#}
# function for mean labels
mean.n <- function(x){
  return(c(y = (median(x))*0.5, label = round(mean(x),2)))
  # experiment with the multiplier to find the perfect position
}

# function for number of observations
give.n <- function(x){
  return(data.frame(y = (median(x)+0.08), label = paste0("n = ", length(x)) ))
  # experiment with the multiplier to find the perfect position
}
#######CN_vs_Biomarker: DEBUG!!, when saved to one pdf, return error?
PlotCN_Biomarker <- function(loop_CNVR, loop_deduped.Intersect.sort.t, loop_deduped.Intersect, One_Biomarker, One_Biomarker_value) {
  # function to plot the box plot for each window with biomaker values
  for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
    # Make sure they merged by rowname?
    one.CN_Biomarker.df <- data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
    one.CN_Biomarker.df$CN <- as.factor(one.CN_Biomarker.df$CN)
    p <- ggplot(one.CN_Biomarker.df, aes_string(x=names(one.CN_Biomarker.df)[1], y=names(one.CN_Biomarker.df)[2])) +
      geom_boxplot(fill = "grey80", colour = "#3366FF") +
      stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
      labs(x = paste0("Window:", loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t), "_winSize:", loop_deduped.Intersect[value,]$Size, "bp" ),
        title = "Boxplot for CN vs associated biomarker:",
        subtitle = paste0("CNVR:", loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "bp", " Biomarker:", as.character(One_Biomarker), "_numberOfWin:", ncol(loop_deduped.Intersect.sort.t) ) )
    #ggsave("Test_PlotCN_Biomarker.png", plot = p)
    print(p)
    }
}

pdf("Test_PlotCN_Biomarker.pdf")
PlotCN_Biomarker(loop_CNVR, loop_deduped.Intersect.sort.t, loop_deduped.Intersect, One_Biomarker, One_Biomarker_value)
dev.off()
########Copy number line plot:
#######CNV_biomarker_rounded_boxplot:
PlotCN_Biomarker_boxplot <- function(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3) {
  BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
  # Sort by sample names' order
  BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
  # select window in the same CNVR
  loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,5] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,6] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,7] == as.numeric(loop_CNVR[3]), ]
  print(paste0("Number of CNVR-window: ",(dim(loop_deduped.Intersect)[1]) ))
  #print(dim(loop_deduped.Intersect))
  # select associate biomarker in the same CNVR
  loop_Biomarker <- unique(Intersect_bed[ Intersect_bed[,7] == as.numeric(loop_CNVR[1]) & Intersect_bed[,8] == as.numeric(loop_CNVR[2]) & Intersect_bed[,9] == as.numeric(loop_CNVR[3]), ][,6])
  for (BiomarkerIndex in 1:length(loop_Biomarker)) {
    One_Biomarker <- loop_Biomarker[BiomarkerIndex]
    One_Biomarker_value <- data.frame(BioMarker.df[, as.character(One_Biomarker), drop = FALSE])
    #colnames(One_Biomarker_value) <- One_Biomarker
    #loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect[,-(1:8)][order(colnames(loop_deduped.Intersect[,-(1:8)]))])
    loop_deduped.Intersect.noCoor <- loop_deduped.Intersect[,-(1:8)]
    loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect.noCoor[,order(colnames(loop_deduped.Intersect.noCoor))])
    # round copy number: if not rounded, will be scattered for not CN=2 regions
    loop_deduped.Intersect.sort.t <- round(loop_deduped.Intersect.sort.t)
    #loop_deduped.Intersect.sort.t <- loop_deduped.Intersect.sort.t[order(rownames(loop_deduped.Intersect.sort.t)),]
    #pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".pdf"))
    pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".rounded.his.boxplot.pdf"))
    PlotCopyNumberSub(loop_CNVR, deduped.Intersect_bed)
    PlotCopyNumber_boxplot(loop_deduped.Intersect)
    #plot CNV_Biomarker boxplot
    for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
      one.CN_Biomarker.df <- data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
      one.CN_Biomarker.df$CN <- as.factor(one.CN_Biomarker.df$CN)
      p <- ggplot(one.CN_Biomarker.df, aes_string(x=names(one.CN_Biomarker.df)[1], y=names(one.CN_Biomarker.df)[2])) +
        geom_boxplot(fill = "grey80", colour = "#3366FF") +
        stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
        labs(x = paste0("CN Window:", loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t), "_winSize:", loop_deduped.Intersect[value,]$Size, "bp" ),
          title = "Boxplot for CN vs associated biomarker:",
          subtitle = paste0("CNVR:", loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "bp", " Biomarker:", as.character(One_Biomarker), "_numberOfWin:", ncol(loop_deduped.Intersect.sort.t) ) )
      #ggsave("myGGplotBoxTest2.png", plot = p)
      print(p)
      }
    dev.off()
  }
}

for ( i in 1:nrow(Sig_CNVR) ) {
  loop_CNVR <- Sig_CNVR[i,]
  PlotCN_Biomarker_boxplot(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3)
}

##debug with 1st CNVR:
loop_CNVR <- Sig_CNVR[1,]
PlotCN_Biomarker_boxplot(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3)


loop_CNVR <- Sig_CNVR[30,]
PlotCN_Biomarker_boxplot(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3)



###try to add summary info into all the CNV_Biomarker boxplot:
summary(one.CN_Biomarker.df)
test <-  summary(one.CN_Biomarker.df)
test[,1]
data.frame(test[,1])


#investigate the merge method for CNV and biomarker tables
data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
merge(data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value["symbol"], by="row.names", all.x=TRUE)
all(rownames(loop_deduped.Intersect.sort.t) == rownames(One_Biomarker_value ))

#ggtitle(paste0("Box plot for CN vs biomarker:", loop_deduped.Intersect[value,1], ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_Size:", loop_deduped.Intersect[value,]$Size, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t) ) ) +
#ggtitle(paste0("Box plot for CN vs biomarker:", loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_Size:", loop_CNVR$CNVRsize, "_Biomarker:", as.character(One_Biomarker), "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t) ) ) +
# plot copy number:

#library("ggplot2")

# the above scripts are for


# error message: > (deduped.Intersect_bed[,c(5:8)] == Sig_CNVR[1,])
#Error in Ops.data.frame(deduped.Intersect_bed[, c(5:8)], Sig_CNVR[1, ]) :
#  ‘==’ only defined for equally-sized data frames

# From the same significant CNVR, make one plot
# ? try: softlink the output to slurm folder for viewing the output in Atom
# Select the pea_3 and the sex+age tables
