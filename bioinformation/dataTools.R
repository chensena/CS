#处理数据的各种小工具集合
#

##按列合并数据
setwd("F:/TRAIN/CHENSEN/genome-photo/chr-link/")

data1 <-  read.table("DfFer-GPX-MT-rename.list", head=T,sep='\t', comment='',quote="")
data2 <-  read.table("DfMT-GPX-Fe-pos.list", head=T,sep='\t', comment='',quote="")
data3 <-  read.table("ChrLen.list", head=T,sep='\t', comment='',quote="")

colnames(data1) <- c("name","altname","model")
colnames(data2) <- c("name","altname")

X <- merge(data2,data1,by = "GID",all=TRUE)
y <- merge(X,data3,by="ChrID")

write.table(y, file = "gene-pos-rename-len.list",sep="\t")

##
#对数转换
setwd("D:\\bash-shell/rna-seq/")
data <-  read.table("D:\\bash-shell/rna-seq/Df_tpm.list", header =T,sep='\t', comment='',quote="")
dim(data)#查看数据维度
#取log值
rownames(data) =data$gene
data=data[,-1] 
logdata <- log(data+1)
write.table(logdata,file= "Df-TPM-log.txt", row.names = T, quote = F, sep="\t")

#排序 筛选

#升序
newdata <- data[order(data$weight),]
newdata

#降序
newdata <- data[order(-data$weight),]
newdata

#按条件取值

chosedata <- data[data$weight >= 0.5,]
dim(chosedata)

write.table(chosedata,file="qu-module10000.txt",row.names = 1,sep="\t")#保存文件

###取数据框前面多少行
chosedata <- newdata[1:10000,]



#突变体数据处理
setwd("E:/XK4-XK5-RNA-SEQ/")

DEG <- read.table("All.DEG_final.xls",header = T,row.names = 1)
colnames(DEG) <- as.character(1:27)

getUncodeGnen <-  function(data){
  list <- vector()#定义一个空的向量
}

#查询STRING 数据库编号

setwd("E:/XK4-XK5-RNA-SEQ/New/DEG-XK5_vs_WT/")
DEGlist <- read.table("DEG.txt",header=T)
protlist <- read.table("E:/STRINGDB/OS-DB/prot-Dfa-39947-last.list",header=T)

#ID查询
merge_list <- merge(DEGlist,protlist,by = "GID")

write.csv(merge_list,file="DEG-STRING-ID.csv")



#相关性图
library(ggplot2)
library(ggpubr)
setwd("F:/TRAIN/CaiXue/")

data <- read.table("CaiXue.txt",header=T)
p <- ggplot(data=data,aes(x=Annual_precipitation,y=MRC)) +
    geom_point(color="red") +
  stat_smooth(method="lm",se=FALSE)+
  stat_cor(method ="pearson")


#过滤低表达样本和基因
library(tidyverse)
library(corrplot)
setwd("F:/TRAIN/CHENSEN/laccase/lignin-fpkm/")
#设置工作目录、
data <- read.table("Filter-DfPOD.tmp.txt",header = T,row.names = 1)

data1 <- data[rowSums(data>=1) >=3,colSums(data>=1)>=3]
write.table(data1,"Filter-POD-FPKM-2.list",sep="\t")


#两组数据的相关性分析
data1  <- read.table("200-400-800-2year.list",header=T)
data2  <- read.table("expersin.list",header=T)
cor_result <- cor(data2,data1,method="spearman")


corrplot(cor_result,addCoef.col="black")
write.table(cor_result,"POD-G-S-H-相关性分析.csv",sep=",")








###############韦恩图绘制
setwd("D:/bio-information/OS/RNA-SEQ/Metal/DEG/")
library(venn)         #韦恩图（venn 包，适用样本数 2-7）
library(VennDiagram) 

# 读取数据文件
venn_dat <- read.delim('DEG-Unit.list')                      # 这里读取了网络上的demo数据，将此处换成你自己电脑里的文件
venn_list <- list(venn_dat[,1], venn_dat[,2], venn_dat[,3], venn_dat[,4], venn_dat[,5], venn_dat[,6], venn_dat[,7],venn_dat[,8])   # 制作韦恩图搜所需要的列表文件
names(venn_list) <- colnames(venn_dat[1:8])    # 把列名赋值给列表的key值

#作图
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1        # 组名字体大小
)

# 更多参数 ?venn查看

# 查看交集详情,并导出结果
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "result.csv", row.names = FALSE, sep = ',', quote = FALSE)




