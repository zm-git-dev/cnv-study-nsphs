# TODO: extend the CNVR regions:
## 26 Aug: extend the window:
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder
#mkdir ExtendWindow
cd /proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow
# CN matrix, whole gnome:
AllCN="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/Select_CNVtable_withName.bed"
# the 39 sifnificant CNVRs
SigCNRV_size="/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/SigCNVR_size_changePvalue.bed"
# extend windows by 2000 bp on each side:?, to do with all regions or just the
cat Chr5_70303300.bed
#5       70303300        70312700        9400
awk '{OFS="\t"; print $1, $2-14000, $3+2600, $4+16600}' Chr5_70303300.bed > Chr5_70303300_17kExtend.bed
#module load bioinfo-tools BEDTools/2.27.1
bedtools intersect -wao -a Chr5_70303300_17kExtend.bed -b $AllCN > AllCN_Chr5_70303300_17kExtend.bed
#less -S AllCN_Chr5_70303300_17kExtend.bed
sort -nk6,6 -nk7,7 AllCN_Chr5_70303300_17kExtend.bed > AllCN_Chr5_70303300_17kExtend.sort.bed
# send to R for making plot
module load bioinfo-tools
module load R/3.4.3
module load R_packages/3.4.3
#in R
list.files()
mybed <- read.table("/proj/sens2016007/nobackup/Zhiwei/SummerProject/FormatCNVbedFolder/estimateP_folder/ExtendWindow/AllCN_Chr5_70303300_17kExtend.sort.bed", header = FALSE, sep="\t")
# remove the laste column, the size of the overlap
# mybedRm <- subset(mybed, select = -c(880))
# pdf("CNVR_Chr5_70303300.pdf")
# #"
# x <- seq(mybed[1,6],mybed[nrow(mybed),7], length = 300)
# y <- seq(-5, 5, length = 300)
# plot(x,y,type="n",xlim=c(mybed[1,6],mybed[nrow(mybed),7]),ylim=c(0,6), xlab = "Chr5:70290700-70315300,24,600bp", ylab = "Copy number")
# for(j in c(1:nrow(mybed))){
#     for (i in 8:(ncol(mybed)-1)) {
#         # RD <- mybed[i]
#         #str(RD)
#         # segments(x0 = t(startPos), y0 = t(RD), x1 = t(endPos), y1 = t(RD), col=col[i])
#         #print(paste0("row and col factor: ", j, " col:", i))
#         #if  (j > 1) {
#           #segments(x0 = mybed[j-1,7], x1 = mybed[j,6],y0 = mybed[j-1,i], y1 = mybed[j,i] )
#         #}
#         segments(x0 = mybed[j,6], x1 = mybed[j,7],y0 = mybed[j,i], y1 = mybed[j,i] )
#     }
# }
# dev.off()
# in ggplot geom_step:
library(ggplot2)
mybed[is.na(mybed)] <- 2
b <- ggplot(mybed)
for (i in 8:(ncol(mybed)-1)) {
  b <- b + geom_step(data=mybed,mapping=aes_string(x=mybed[,6], y=mybed[,i]))
}
ggsave("CNVR_Chr5_70303300_ggplot.pdf")

b <- ggplot(mybed)
for (i in 8:30) {
  b + geom_step(data=mybed,mapping=aes_string(x=mybed[,6], y=mybed[,i]))
}
ggsave("CNVR_Chr5_70303300_ggplot.pdf")

pdf("CNVR_Chr5_70303300_ggplot.pdf") #"
print(b)
dev.off()
b <- ggplot(mybed) + geom_step(data=mybed, mapping=aes_string(x=mybed[,6], y=mybed[,i]))
b + geom_step(data=mybed, mapping=aes_string(x=mybed[,6], y=mybed[,i+1]))

library(tidyr)
#remove the intersect size info:
mybedRm <- subset(mybed, select = -c(880))
long_mybed <- mybedRm %>% gather(Sample, CN, 8:ncol(mybedRm))

#ggplot(long_mybed)

b <- ggplot(long_mybed,aes(V6,CN))
b + geom_step(aes(group = Sample))
ggsave("CNVR_Chr5_70303300_ggplot.pdf") #"
##GGplot step function
b <- ggplot(long_mybed)
b <- b + geom_step(aes(x=V6,y=CN, group=Sample ))
b + labs(x = "Position", y = "Copy Number", title = "Step Plot")
ggsave("CNVR_Chr5_70303300_ggplot.pdf")
