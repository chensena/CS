#TPM值的计算，FPKM值的计算

rm(list = ls())
options(stringsAsFactors = F)
rawcount <- read.csv("C:/Users/Desktop/2cell.id.txt",skip = 1,sep="\t",row.names = 1)
rawcount$gene_name <- rownames(rawcount)
length <- read.csv('C:/Users/Desktop/exons_gene_lens_df_ensembl（小鼠）(1).csv ', header = T, row.names = 1)
colnames(length) <- c('gene_name', 'length')

#install.packages("dplyr")
library(dplyr)
length <- distinct(length, gene_name, .keep_all = T)
mergecount <- merge(rawcount, length, by = 'gene_name')
FPKMlength <- mergecount[,c(1, ncol(mergecount))]
FPKMcount <- mergecount[, -c(2, ncol(mergecount))]
rownames(FPKMcount) <- FPKMcount$gene_name
FPKMcount <- FPKMcount[,-1]

## 计算 TPM ##

kb <- FPKMcount$Length / 1000
kb
countdata <- FPKMcount[,5:7]   #r的索引是从1开始的，5：7选择的是count里面每个样本对应的reads数的列
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
#将上面计算好的tpm保存到本地
ootpm <- as.data.frame(tpm)
write.csv(ootpm, file="C:/Users/Desktop/ootpm.csv",quote=FALSE)


## 计算 FPKM ##
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
ofpkm <- as.data.frame(fpkm)
write.csv(ofpkm , file="C:/Users/Desktop/2tpm.csv",quote=FALSE)
