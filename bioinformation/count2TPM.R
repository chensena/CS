library(dplyr)
setwd("D:\\bio-information/MA/rna_seq/genomic/")
length=read.table("MA_length.txt",header=T,sep='\t',check.names=F)
#提取基因长度
count <- read.table(file="D:/bio-information/MA/rna_seq/SRR/out/SAM/Dla_All_Count.txt",header = T,sep = "\t")

merge<-left_join(count,length,by="Gene")#根据基因那列进行合并
merge <- na.omit(merge)#删除错误值行
write.csv(merge,file = "merge.csv",sep = "\t")#读出文件，直接往下运行也许

mycounts<-read.csv("merge.csv")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)#最后一列Length是基因长度

#TPM计算
kb <- mycounts$Length / 1000
kb
countdata <- mycounts[,1:24]
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.table(tpm,file="tpm.xls",sep="\t",quote=F)