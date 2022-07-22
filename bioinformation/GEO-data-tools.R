#GEO数据库挖掘

#数据下载
###########################################
# GEO accession : GSE42872
# Platforms     : GPL6244
# BioProject    : PRJNA183688
##########################################

setwd("D:/bio-information/OS/RNA-SEQ/ALL-RNA-SEQ/GSE40549/")

#####数据下载#####
if(!require(GEOquery)) BiocManager::install("GEOquery") # 安装包
package.version("GEOquery") # 查看版本
help(package = "GEOquery") # 查看GEOquery中的函数
library(GEOquery) # 加载包
library(tidyverse)
library(genefilter)

exprset <- read.table(file = 'GSE40549_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      check.names = F) #读取表达数据

#循环去除矩阵引号
## nrow(AA)-----这个矩阵的行数
for (i in 1:nrow(exprset) ){
  x=exprset[i,1]  # 赋值
  x=as.character(x) #化作字符串
  a=gsub('["]', '', x)  #去双引号
  exprset[i,1]=a  #给矩阵重新赋值
}

#去除第一行的引号

for (i in 1:ncol(exprset)){
  x=exprset[1,i]  # 赋值
  x=as.character(x) #化作字符串
  a=gsub('["]', '', x)  #去双引号
  exprset[1,i]=a  #给矩阵重新赋值
}
#定义列名
colnames(exprset) = exprset[1,]

exprset2 <- exprset[-1,]  #
去除第一列数据
#####
#芯片ID转换


id_table <-read.table("GPL6864_old_no_an.txt",header=TRUE,sep="\t",na.strings = "",quote="")
colnames(id_table)
require(data.table)
probe2symbol <- id_table[,c("ID","Accessions","GB_ACC")]

#得到ID转换矩阵


###矩阵合并
merged_expr_df <- merge(x=exprset2,y=probe2symbol,by.x="ID_REF",by.y="ID")

#去空白值

#merged_expr_df[merged_expr_df == ""] <- NA
#空白值用NA替代，去除Accessions列为NA的行
filt_merged_expr_df <- merged_expr_df[-which(is.na(merged_expr_df$Accessions)),]
#删除不需要的行
next_exprset <- select(filt_merged_expr_df,-c("ID_REF","GB_ACC"))
#分列，Accessions分列为“symbol”,"a","b","c"
split_last_exprset <- separate(data = last_exprset, col = Accessions, into = c("symbol", "a","b","c"), sep = "\\|")
last_exprset <- select(split_last_exprset,-c("a","b","c"))
#得到表达矩阵
#去除重复探针,重复值取最大值


uniq_exprset <-  aggregate(x = last_exprset,by = list(last_exprset$symbol), FUN = max)

row.names(uniq_exprset) <- uniq_exprset$symbol#symbol列变为列名

output_exprset <- select(uniq_exprset,-c("Group.1","symbol"))#去除不需要的列

write.csv(output_exprset,"GSE40549-last-exprset.csv")








####表达矩阵批量下载器

library(GEOquery)
library(tidyverse)
sewd("D:/bio-information/OS/RNA-SEQ/Metal/")

GSE_list <- c("GSE33375","GSE33376","GSE34895","GSE35502","GSE36760","GSE41719","GSE41733","GSE63152")

for (i in GSE_list){
  GSE_id <- i
  gset <- getGEO(GSE_id,destdir = ".",
                 AnnotGPL = F,getGPL = F)
  exprset <- data.frame(exprs(gset[[1]])) # 推荐
  write.csv(exprset,paste0(GSE_id,".csv"))
}

#网络不通畅--手动下载
#GSE63152
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE63152/")

exprset <- read.table(file = 'GSE63152_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据


#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))



#取每行log2FC的平均值
last_exprset$row_mean  <- apply(last_exprset,1,mean)
last_exprset$change = as.factor(ifelse(abs(last_exprset$row_mean) > 1.5,ifelse(last_exprset$row_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE63152-exprset.csv")



#GSE41733
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE41733/")


exprset <- read.table(file = 'GSE41733_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据

#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))


last_exprset$short_mean <- apply(last_exprset[,1:3],1,mean)
last_exprset$long_mean <- apply(last_exprset[,4:6],1,mean)
last_exprset$S_change = as.factor(ifelse(abs(last_exprset$short_mean) > 1.5,ifelse(last_exprset$short_mean > 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$L_change = as.factor(ifelse(abs(last_exprset$long_mean) > 1.5,ifelse(last_exprset$long_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE41733_exprset.csv")



#GSE41719
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE41719/")
exprset <- read.table(file = 'GSE41719_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据

#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))


last_exprset$short_mean <- apply(last_exprset[,1:3],1,mean)
last_exprset$long_mean <- apply(last_exprset[,4:6],1,mean)
last_exprset$S_change = as.factor(ifelse(abs(last_exprset$short_mean) > 1.5,ifelse(last_exprset$short_mean > 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$L_change = as.factor(ifelse(abs(last_exprset$long_mean) > 1.5,ifelse(last_exprset$long_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE41719_exprset.csv")



#GSE36760
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE36760/")


exprset <- read.table(file = 'GSE36760_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据


#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))



#取每行log2FC的平均值
last_exprset$row_mean  <- apply(last_exprset,1,mean)
last_exprset$change = as.factor(ifelse(abs(last_exprset$row_mean) > 1.5,ifelse(last_exprset$row_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE36760_exprset.csv")







#GSE35502
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE35502/")


exprset <- read.table(file = 'GSE35502_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据


#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))


#分析上调和下调基因
last_exprset$iR_WR_mean <- apply(last_exprset[,1:6],1,mean)
last_exprset$KS_WS_mean <- apply(last_exprset[,7:12],1,mean)
last_exprset$WR_mean <- apply(last_exprset[,13:16],1,mean)
last_exprset$Ws_mean <- apply(last_exprset[,17:20],1,mean)

last_exprset$iR_WR_change = as.factor(ifelse(abs(last_exprset$iR_WR_mean) > 1.5,ifelse(last_exprset$iR_WR_mean> 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$KS_WS_change = as.factor(ifelse(abs(last_exprset$KS_WS_mean) > 1.5,ifelse(last_exprset$KS_WS_mean > 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$WR_change = as.factor(ifelse(abs(last_exprset$WR_mean) > 1.5,ifelse(last_exprset$WR_mean > 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$WS_change = as.factor(ifelse(abs(last_exprset$Ws_mean) > 1.5,ifelse(last_exprset$Ws_mean > 1.5 ,'UP','DOWN'),'STABLE'))


write.csv(last_exprset,"GSE35502_exprset.csv")

#CSE34895
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE34895/")

exprset <- read.table(file = 'GSE34895_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据

#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))


last_exprset$Cu_mean <- apply(last_exprset[,1:3],1,mean)
last_exprset$Cd_mean <- apply(last_exprset[,4:6],1,mean)
last_exprset$Cu_change = as.factor(ifelse(abs(last_exprset$Cu_mean) > 1.5,ifelse(last_exprset$Cu_mean > 1.5 ,'UP','DOWN'),'STABLE'))
last_exprset$Cd_change = as.factor(ifelse(abs(last_exprset$Cd_mean) > 1.5,ifelse(last_exprset$Cd_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE34895_exprset.csv")




#GSE33376
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE33376/")

exprset <- read.table(file = 'GSE33376_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据


#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))



#取每行log2FC的平均值
last_exprset$row_mean  <- apply(last_exprset,1,mean)
last_exprset$change = as.factor(ifelse(abs(last_exprset$row_mean) > 1.5,ifelse(last_exprset$row_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE33376-exprset.csv")



#GSE33375
setwd("D:/bio-information/OS/RNA-SEQ/Metal/GSE33375/")

exprset <- read.table(file = 'GSE33375_series_matrix.txt',
                      sep = '\t',
                      quote = '',
                      fill = T, 
                      comment.char = "!",
                      header=T,
                      check.names = F) #读取表达数据


#数据分列
split_exprset <- separate(data = exprset, col = ID_REF, into = c("symbol", "a","b","c"), sep = "\\|")
#重复值取最大值

symbol_exprset <- select(split_exprset,-c("a","b","c"))
#获得有重复值的表达矩阵
uniq_exprset <-  aggregate(x = symbol_exprset,by = list(symbol_exprset$symbol), FUN = mean)

row.names(uniq_exprset) <- uniq_exprset$Group.1

last_exprset <- select(uniq_exprset,-c("symbol","Group.1"))



#取每行log2FC的平均值
last_exprset$row_mean  <- apply(last_exprset,1,mean)
last_exprset$change = as.factor(ifelse(abs(last_exprset$row_mean) > 1.5,ifelse(last_exprset$row_mean > 1.5 ,'UP','DOWN'),'STABLE'))
write.csv(last_exprset,"GSE33375-exprset.csv")






