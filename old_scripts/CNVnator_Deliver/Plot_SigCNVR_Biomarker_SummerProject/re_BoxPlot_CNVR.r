### reorganize the script to debug the plot function in SignificantCNVR.r
# Files directory and loaded software:
# in Bash
#cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3

# in R:
#rm(list=ls())
library(ggplot2)

# in Bash
#less -S SigCNVR_SigAssociation_CN.bed

#load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
# ls()
#AgeSex <- read.csv(file="/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/Physical/KA06_KA09_age_sex.csv", header=TRUE, sep=",")
#Intersect_bed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/SigCNVR_SigAssociation_CN.bed", header = FALSE, check.names=FALSE, sep = "\t")
# !remove the last 2 columns from bedtools intersect overlap region SIZE, add header to the data frame.
#Intersect_bed <- subset(Intersect_bed, select = -c(887,888))
# remove redudant window coordiante, proof by:
#> a1 <- data.frame(Intersect_bed[,c(5,6,7)])
#> a2 <- data.frame(Intersect_bed[,c(12,13,14)])
#> colnames(a1) <- c("1", "2", "3")
#> colnames(a2) <- c("1", "2", "3")
#> identical(a1,a2)
#[1] TRUE
#remove the extra window coordinate mentioned above:
#Intersect_bed <- subset(Intersect_bed, select = -c(12,13,14))

#CNVbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.txt", sep="\t", header = TRUE, comment.char = "", check.names = FALSE)
#SampleName <- colnames(CNVbed)
#rm(CNVbed)

#c("CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[1:3], "WinSize", "P_value", "StatsInfo", "Biomarker", SampleName[4:length(SampleName)])
# add column names to the intersect bed:
#colnames(Intersect_bed) <- c("CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[1:3], "WinSize", "P_value", "StatsInfo", "Biomarker", SampleName[4:length(SampleName)])
#Intersect_bed[1:5,1:20]
#BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
# Sort by sample names' order
#BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
#summary(colnames(Intersect_bed)[12:ncol(Intersect_bed)]==rownames(BioMarker.df))

#########function:
# Select CNVR:
SelectCNVR_function <- function(Intersect_bed) {
  # select the 1~4 columns: coordinate and Size of CNVRs, use unique function to have the CNVRs appear once:
  # e.g: Sig_CNVR <- SelectCNVR_function(Intersect_bed)
  Sig_CNVR <- unique (Intersect_bed[,1:4])
  return(Sig_CNVR)
}
#Design: loop to generate the figure for the 30 regions:
#for ( i in 1:nrow(Sig_CNVR) ) {
#  loop_CNVR <- Sig_CNVR[i,]
#  #PlotCN_Biomarker_boxplot(loop_CNVR, deduped.Intersect_bed, Intersect_bed, pea_3)
#}
# Generate dedupe intersect_bed, (no duplicate by multiple association to biomarkers)
GenerateDeduptIntersect_bed <- function(Intersect_bed) {
  # generate the unique list of CNVR-windw coordinate + CN info for all sample:
  #e.g deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)
  deduped.Intersect_bed <- unique(Intersect_bed[, c(1:8,12:ncol(Intersect_bed))])
  return(deduped.Intersect_bed)
}
# Generate regions in the same CNVR, only corr and CN values:
GenerateLoopDedupedIntersect_bed <- function(loop_CNVR, deduped.Intersect_bed) {
  # given one CNVR at a time and get the windows within the CNVR iwht CN info
  #e.g:loop_deduped.Intersect <- GenerateLoopDedupedIntersect_bed(loop_CNVR, deduped.Intersect_bed)
  loop_deduped.Intersect <- deduped.Intersect_bed[deduped.Intersect_bed[,"CNVRchr"] == as.numeric(loop_CNVR[1]) & deduped.Intersect_bed[,"CNVRstart"] == as.numeric(loop_CNVR[2]) & deduped.Intersect_bed[,"CNVRend"] == as.numeric(loop_CNVR[3]), ]
  return(loop_deduped.Intersect)
}
# Given loop_CNVR: select Biomarker(s) associated with the same CNVR
SelectBiomarker_function <- function(loop_CNVR, Intersect_bed) {
  # output the associated biomarker(s) in the given CNVR
  #e.g: loop_Biomarker <- SelectBiomarker_function(loop_CNVR, Intersect_bed)
  loop_Biomarker <- unique(Intersect_bed[ Intersect_bed[,"CNVRchr"] == as.numeric(loop_CNVR[1]) & Intersect_bed[,"CNVRstart"] == as.numeric(loop_CNVR[2]) & Intersect_bed[,"CNVRend"] == as.numeric(loop_CNVR[3]), ][,"Biomarker"])
  return(loop_Biomarker)
}
# then: for (BiomarkerIndex in 1:length(loop_Biomarker)) {...}
# Generate One_Biomarker (name) and One_Biomarker_value (population value)
# BiomarkerIndex <- 1
Generate_OneBiomarkerValue <- function (BioMarker.df, loop_Biomarker, BiomarkerIndex) {
  # generete data for one biomarker value given the BiomarkerIndex in loop_Biomarker:
  # e.g One_Biomarker_value <- Generate_OneBiomarkerValue(BioMarker.df, loop_Biomarker, BiomarkerIndex)
  One_Biomarker <- loop_Biomarker[BiomarkerIndex]
  One_Biomarker_value <- data.frame(BioMarker.df[, as.character(One_Biomarker), drop = FALSE]) #!!!error!# no column name? # if not as.character, it would use the index(numeric!) insize the vector
  #colnames(One_Biomarker_value) <- One_Biomarker
  return(One_Biomarker_value)
}

Generate_transpose_noCoorCN <- function (loop_deduped.Intersect) {
  # remove the coordinate of the CNVR and window in the CN table
  # sort by sample name?
  # e.g loop_deduped.Intersect.sort.t <- Generate_transpose_noCoorCN(loop_deduped.Intersect)
  loop_deduped.Intersect.noCoor <- loop_deduped.Intersect[,-(1:8)]
  loop_deduped.Intersect.sort.t <- t(loop_deduped.Intersect.noCoor[,order(colnames(loop_deduped.Intersect.noCoor))])
  # round copy number: if not rounded, will be scattered for not CN=2 regions
  loop_deduped.Intersect.sort.t <- round(loop_deduped.Intersect.sort.t)
  loop_deduped.Intersect.sort.t <- loop_deduped.Intersect.sort.t[order(rownames(loop_deduped.Intersect.sort.t)),]
  # format as data frame, otherwise in cases of only one window, return error on ncol(loop_deduped.Intersect.sort.t)
  return(data.frame(loop_deduped.Intersect.sort.t))
}

#loop_CNVR <- Sig_CNVR[1,]
#deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)
#loop_deduped.Intersect <- GenerateLoopDedupedIntersect_bed(loop_CNVR, deduped.Intersect_bed)
#loop_Biomarker <- SelectBiomarker_function(loop_CNVR, Intersect_bed)
#BiomarkerIndex <- 1
#One_Biomarker_value <- Generate_OneBiomarkerValue(BioMarker.df, loop_Biomarker, BiomarkerIndex)
#loop_deduped.Intersect.sort.t <- Generate_transpose_noCoorCN(loop_deduped.Intersect)

##test:
#> loop_CNVR
#  CNVRchr CNRVStart   CNVRend CNVRsize
#1       1 161640580 161642980     2400
#> dim(deduped.Intersect_bed)
#[1] 285 880
#> dim(loop_deduped.Intersect)
#[1]   7 880
#> loop_Biomarker
#[1] ONC2_194_FCRLB
#17 Levels: CVD2_133_IL-18 CVD2_184_PD-L2 CVD2_186_hOSCAR ... ONC2_194_FCRLB
#> dim(One_Biomarker_value)
#[1] 872   1
#> dim(loop_deduped.Intersect.sort.t)
#[1] 872   7

### plot function:
# plot histogram for each CNVR, (merge info of all windows in the CNVR)
PlotCopyNumberHis <- function(loop_deduped.Intersect) {
  # Extract all CN values in all windows:
  CopyNumber <- c(t(loop_deduped.Intersect[,9:ncol(loop_deduped.Intersect)]))
  CN.df <- data.frame(CopyNumber)
  CN.p <- ggplot(data= CN.df , aes(x= CopyNumber)) +
     geom_histogram() +
     ggtitle(paste0("CN histogram for CNVR:", loop_deduped.Intersect$CNVRchr[1], ":", loop_deduped.Intersect$CNVRstart[1], "-", loop_deduped.Intersect$CNVRend[1], "_Size:", loop_deduped.Intersect$CNVRsize[1], "bp", "_NoWin:", nrow(loop_deduped.Intersect) )) +
     labs(x = "Copy Number")
  print(CN.p)
}
# test:
#library(ggplot2)
#pdf("testPlotCNhisto.pdf")
#PlotCopyNumberHis(loop_deduped.Intersect)
#dev.off()

# plot boxplot for each CNVR, one box refer to one window in the CNVR:
PlotCopyNumber_boxplot <- function (loop_deduped.Intersect) {
  Cor.CV.table = data.frame()
  for (rowIndex in 1:nrow(loop_deduped.Intersect)) {
    # only pick the starting position of a window
    Cor.CV.df <- data.frame( rep(loop_deduped.Intersect[rowIndex,"Star"], ncol(loop_deduped.Intersect)-8 ), t(loop_deduped.Intersect[rowIndex, 9:ncol(loop_deduped.Intersect)]) )
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
# test
#pdf("testPlotCNboxplot.pdf")
#PlotCopyNumber_boxplot(loop_deduped.Intersect)
#dev.off()

# plot boxplot for window vs biomarker values, group as per CNVR:
# function for mean labels
mean.n <- function(x){
  return(data.frame(y = (median(x)-0.09), label = paste0("mean= ", round(mean(x),2)) ) )
  # experiment with the multiplier to find the perfect position
}
# function for number of observations
give.n <- function(x){
  return(data.frame(y = (median(x)+0.09), label = paste0("n = ", length(x)) ))
  # experiment with the multiplier to find the perfect position
}

PlotCN_Biomarker <- function(loop_deduped.Intersect.sort.t, loop_deduped.Intersect, One_Biomarker, One_Biomarker_value) {
  # function to plot the box plot for each window with biomaker values
  print(paste0("Number of CNVR-window: ",(dim(loop_deduped.Intersect.sort.t)[2]) ))
  for (value in 1:ncol(loop_deduped.Intersect.sort.t)) {
    # Make sure they merged by rowname?
    one.CN_Biomarker.df <- data.frame( data.frame(CN = round(loop_deduped.Intersect.sort.t[,value])), One_Biomarker_value)
    one.CN_Biomarker.df$CN <- as.factor(one.CN_Biomarker.df$CN)
    # remove individuals that have NA copy number:
    one.CN_Biomarker.df <- one.CN_Biomarker.df[!is.na(one.CN_Biomarker.df$CN),]
    # remove rows with NA in CN column
    p <- ggplot(one.CN_Biomarker.df, aes_string(x=names(one.CN_Biomarker.df)[1], y=names(one.CN_Biomarker.df)[2])) +
      geom_boxplot(fill = "grey80", colour = "#3366FF") +
      stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
      stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red") +
      labs(x = paste0("Copy Number in Window:", loop_deduped.Intersect[value,]$Chr, ":", loop_deduped.Intersect[value,]$Star, "-", loop_deduped.Intersect[value,]$End, "_No.", value, "/", ncol(loop_deduped.Intersect.sort.t), "_winSize:", loop_deduped.Intersect[value,]$WinSize, "bp" ),
        title = "Boxplot for CN vs associated biomarker:",
        subtitle = paste0("CNVR:", loop_deduped.Intersect$CNVRchr[1], ":", loop_deduped.Intersect$CNVRstart[1], "-", loop_deduped.Intersect$CNVRend[1], "_Size:", loop_deduped.Intersect$CNVRsize[1], "bp", " Biomarker:", as.character(One_Biomarker), "_numberOfWin:", ncol(loop_deduped.Intersect.sort.t) ) )
    #ggsave("Test_PlotCN_Biomarker.png", plot = p)
    print(p)
    }
}

#pdf("testPlotCN_Biomarker_boxplot.pdf")
#PlotCN_Biomarker(loop_deduped.Intersect.sort.t, loop_deduped.Intersect, One_Biomarker, One_Biomarker_value)
#dev.off()

# build function to show all:
#loop_CNVR <- Sig_CNVR[1,]
#deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)
#loop_deduped.Intersect <- GenerateLoopDedupedIntersect_bed(loop_CNVR, deduped.Intersect_bed)
#loop_Biomarker <- SelectBiomarker_function(loop_CNVR, Intersect_bed)
#BiomarkerIndex <- 1
#One_Biomarker_value <- Generate_OneBiomarkerValue(BioMarker.df, loop_Biomarker, BiomarkerIndex)
#loop_deduped.Intersect.sort.t <- Generate_transpose_noCoorCN(loop_deduped.Intersect)

######main function:
load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
Intersect_bed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/SigCNVR_SigAssociation_CN.bed", header = FALSE, check.names=FALSE, sep = "\t")
# !remove the last 2 columns from bedtools intersect overlap region SIZE, add header to the data frame.
Intersect_bed <- subset(Intersect_bed, select = -c(887,888))
#remove the extra window coordinate mentioned above:
Intersect_bed <- subset(Intersect_bed, select = -c(12,13,14))
# read sample names:
CNVbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.txt", sep="\t", header = TRUE, comment.char = "", check.names = FALSE)
SampleName <- colnames(CNVbed)
rm(CNVbed)
# give name to the Intersect_bed
# add column names to the intersect bed:
colnames(Intersect_bed) <- c("CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[1:3], "WinSize", "P_value", "StatsInfo", "Biomarker", SampleName[4:length(SampleName)])
# select biomarker:
BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
# Sort by sample names' order
BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)
Sig_CNVR <- SelectCNVR_function(Intersect_bed)
# check sample name in correct order:
#summary(colnames(Intersect_bed)[12:ncol(Intersect_bed)]==rownames(BioMarker.df))
####### Plot main function:
#Sig_CNVR <- SelectCNVR_function(Intersect_bed)
#loop_CNVR <- Sig_CNVR[1,] # test with the first CNVR
#deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)

PlotMain <- function(loop_CNVR, deduped.Intersect_bed, Intersect_bed, BioMarker.df) {
  loop_deduped.Intersect <- GenerateLoopDedupedIntersect_bed(loop_CNVR, deduped.Intersect_bed)
  loop_Biomarker <- SelectBiomarker_function(loop_CNVR, Intersect_bed)
  for (BiomarkerIndex in 1:length(loop_Biomarker)) {
    One_Biomarker <- loop_Biomarker[BiomarkerIndex]
    One_Biomarker_value <- Generate_OneBiomarkerValue(BioMarker.df, loop_Biomarker, BiomarkerIndex)
    loop_deduped.Intersect.sort.t <- Generate_transpose_noCoorCN(loop_deduped.Intersect)
    # the plots will be saved in the same directory:
    pdf(paste0(loop_CNVR$CNVRchr, ":", loop_CNVR$CNVRstart, "-", loop_CNVR$CNVRend, "_", as.character(One_Biomarker), ".rounded.his.boxplot.pdf"))
    PlotCopyNumberHis(loop_deduped.Intersect)
    PlotCopyNumber_boxplot(loop_deduped.Intersect)
    PlotCN_Biomarker(loop_deduped.Intersect.sort.t, loop_deduped.Intersect, One_Biomarker, One_Biomarker_value)
    dev.off()
  }
}
# test:
#loop_CNVR <- Sig_CNVR[1,]
#PlotMain(loop_CNVR, deduped.Intersect_bed, Intersect_bed, BioMarker.df)

# Main run
for (i in 1:nrow(Sig_CNVR)) {
  loop_CNVR <- Sig_CNVR[i,]
  PlotMain(loop_CNVR, deduped.Intersect_bed, Intersect_bed, BioMarker.df)
}

##29 Aug, change P value:
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/Figures_changeP
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3

load("/proj/sens2016007/nobackup/NSPHS_phenotype_data/Phenotypes/PEA/pea3.rntransformed.RData")
Intersect_bed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_SigAssociation_CN_changeP.bed", header = FALSE, check.names=FALSE, sep = "\t")
# !remove the last 2 columns from bedtools intersect overlap region SIZE, add header to the data frame.
Intersect_bed <- subset(Intersect_bed, select = -c(887,888))
#remove the extra window coordinate mentioned above:
Intersect_bed <- subset(Intersect_bed, select = -c(12,13,14))
# read sample names:
CNVbed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.txt", sep="\t", header = TRUE, comment.char = "", check.names = FALSE)
SampleName <- colnames(CNVbed)
rm(CNVbed)
# give name to the Intersect_bed
# add column names to the intersect bed:
colnames(Intersect_bed) <- c("CNVRchr", "CNVRstart", "CNVRend", "CNVRsize", SampleName[1:3], "WinSize", "P_value", "StatsInfo", "Biomarker", SampleName[4:length(SampleName)])
# select biomarker:
BioMarker.df <- pea_3[rownames(pea_3) %in% colnames(Intersect_bed), ]
# Sort by sample names' order
BioMarker.df <- BioMarker.df[order(rownames(BioMarker.df)), ]
deduped.Intersect_bed <- GenerateDeduptIntersect_bed(Intersect_bed)
Sig_CNVR <- SelectCNVR_function(Intersect_bed)
for (i in 1:nrow(Sig_CNVR)) {
  loop_CNVR <- Sig_CNVR[i,]
  PlotMain(loop_CNVR, deduped.Intersect_bed, Intersect_bed, BioMarker.df)
}
