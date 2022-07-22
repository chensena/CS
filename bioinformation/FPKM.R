#!/usr/bin/R
setwd("")
library(GenomicFeatures)
setwd("D:\\bio-information\\RNA-SEQ\\FPKM-CAU")
txdb <- makeTxDbFromGFF("D:\\bio-information\\compare-rna-seq\\DF-genenomic\\DFA.gff",format="gff")
#创建txDB对象，输入文件可以是GTF,GFF3
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
n=t(as.data.frame(exons_gene_lens))
write.table(n,'DFA_length.txt',col.names=F,row.names=T,quote=F,sep='\t')
length=read.table("Maize_length.txt",header=F,sep='\t',check.names=F)
names(length)<-c('gene','length')

#导入所有reads的count值,对特殊符号进行处理
count <- read.table(file="Maize-count.txt",header = T,sep = "\t")

count[,1] = gsub(':','.',count[,1])
count[,1] = gsub('-','.',count[,1])
count[,1] = gsub(' ','.',count[,1])
count[,1] = gsub('/','.',count[,1])

merge<-merge(count,length,by = 'gene') #匹配gene_count与gene_length
dim(merge)
count <- merge[,1:(dim(merge)[2]-1)]
gene_num <- dim(merge)[1]
sample_num <- dim(merge)[2]-2 #减去gene_name和gene_length列
#从第二列开始，将每个样本的count循环转化成FPKM
i <- 2

repeat{
  mapped_reads <- sum(merge[1:gene_num,i])#计算每个样本的mapped reads数
  FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])#计算FPKM值
  FPKM <- pmax(FPKM,0)#去掉矫正带来的负值
  count = data.frame(count[1:gene_num,],FPKM)#添加到count表格后面i
  i <- i + 1
  
  if(i > sample_num+1){
    break
  }
}

#生成表格列名称
head(count)
count_colname <- read.table("Maize-count.txt",header = T,nrow = 1,as.is=TRUE)
FPKM_colname <- paste(count_colname[1,],"_FPKM",sep="")
FPKM_colname
colname <- c(count_colname,FPKM_colname)
col_name <- colname[-which(colname=="gene_FPKM")]#删除gene_FPKM
col_name
names(count) <- col_name
head(count)

#生成表格
write.csv(count,"Maize_FPKM.csv")
write.table(count[,c(1,(sample_num+2):(sample_num*2+1))],"Maize_FPKM_table.txt",row.names = FALSE, quote = FALSE, sep = "\t")
#/usr/bin/bash sed 's/,/\t/g' FPKM_table_tmp.txt > FPKM_table.txt
head(read.csv("Maize_FPKM.csv"))