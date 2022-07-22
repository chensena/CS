#tpm 做 TMM均一化
library(edgeR)
library(ggplot2)
library(pheatmap)

setwd("D:\\bio-information/compare-rna-seq/OS-ZM-DF-DFA-DFB-DFC/OS-ZM-DFA/")
data<-read.delim("D:\\bio-information/compare-rna-seq/OS-ZM-DF-DFA-DFB-DFC/OS-ZM-DFA/OS-ZM-A.list",row.names = 1,sep='\t',check.names=FALSE)
group <- factor(rep(c('OS','ZM','DFA'), each = 3))
data <- as.matrix(data)

setwd("D:\\bio-information\\compare-rna-seq\\ALL-FPKM\\one-to-one TPM")
data <- read.delim('D:\\bio-information\\compare-rna-seq\\ALL-FPKM\\one-to-one TPM\\all-121-tpm-notitle.txt', row.names = 1, sep = '\t', check.names = FALSE)
group <- factor(rep(c('C10', 'C100', 'C200','C400','C50','C800'), each = 3))
group <- factor(rep(c('Branch','C10','C100','C200','C400','C50','C800','CSh','Inod','Labud','MLeaf','Node','Rhbud','Rhne','Root','ShB','Yleaf'),each = 3))
#
group <- factor(colnames(data))
group<-factor(data)



setwd("D:/bio-information/comepare-genomic/OS-MAIZE-DFA-MAA/OS-ZM-DFA-MAA-121/")
data<-read.delim("D:\\bio-information/comepare-genomic/OS-MAIZE-DFA-MAA/OS-ZM-DFA-MAA-121/OS-ZM-DFA-MAA-leaf-tmp.txt",row.names = 1,sep='\t',check.names=FALSE)
group <- factor(rep(c('OS','DFA','MAA','ZM'), each = 3))


#不同阶段的笋子的均一化
setwd("D:/bio-information/WGCNA/SUN/")
data<-read.delim("D:/bio-information/WGCNA/SUN/Df_SunZi_tpm.list",row.names = 1,sep='\t',check.names=FALSE)
group <- factor(rep(c('C10','C100','C200','C400','C50','C800'), each = 3))


#梁山慈竹-麻竹均一化

setwd("D:/bio-information/comepare-genomic/Dfa-MA-orthofinder/ABC/DEG/")
data<-read.delim("MA-DF-Tmp-121.txt",row.names = 1,sep='\t',check.names=FALSE)
data <- data[rowMeans(data) >1,]
sample <- read.table("MA-DF-simple.list",header=T,row.names = 1,sep="\t")
group <- as.factor(sample$type)


#麻竹均一化
setwd("D:\\bio-information\\WGCNA/MA-WGCNA/TPM/")
data<-read.delim("Ma_all_DEG_tpm.csv",row.names = 1,sep=',',check.names=FALSE)
data <- data[rowMeans(count) >1,]
sample <- read.table("D:/bio-information/MA/rna_seq/DEseq2/sample.list",header=T,row.names = 1,sep="\t")
group <- as.factor(sample$type)


#麻竹所有组织均一化
setwd("D:/bio-information/comepare-genomic/Dfa-MA-orthofinder/ABC/")
data<-read.delim("DF-MA-tpm-tmm.csv",row.names=1,sep=',',check.names=FALSE)
group <- factor(rep(c('Branch','C10','C100','C200','C400','C50','C800','CSh','Inod','Labud','MLeaf','Node','Rhbud','Rhne','Root','ShB','Yleaf','abshoot','internode','lateral_bud','leaf','root','shoot','shoot_shell','stem'),each = 3))
group <- as.factor(group)



data <- as.matrix(data)
data <- as.numeric(data)
y <- DGEList(counts = data, group = group)
y <- calcNormFactors(y)
logcpm <- cpm(y, prior.count=3, log=TRUE)
write.table(logcpm, file="all_leaf_TMM.xls", sep="\t", quote=F, row.names=T, col.names=T)



dgList <- estimateCommonDisp(y,verbose=TRUE)
dgList <- estimateTagwiseDisp(dgList,verbose=TRUE)
norm_counts.table <- t(t(dgList$pseudo.counts)*(dgList$samples$norm.factors))
write.table(norm_counts.table, file="DF-MA-_tpm-tmm-tmm.txt", sep="\t", quote=F)
pheatmap(log(data+1),cluster_rows=T,cluster_cols=T,scale="none",border_color="white",color=colorRampPalette(rev(c("red","white","blue")))(102))







####求平均表达量
setwd("D:/bash-shell/rna-seq/")

myTpm <- read.table("DF_tpm-tmm.txt",header = T,row.names = 1)



group <- factor(rep(c('Branch','C10','C100','C200','C400','C50','C800','CSh','Inod','Labud','MLeaf','Node','Rhbud','Rhne','Root','ShB','Yleaf'),each = 3))


myMeanFun<-function(x){
  tapply(as.double(x),group,mean)
  
}

meanTpm <- t(apply(myTpm,1,myMeanFun))  #计算三列的均值

write.csv(meanTpm,"DF-meanTpm-tmm.csv")



#求平均表达量，mini版

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(export)
library(png)


setwd("E:/Sub_chr_E/expersion/ABC-expersion/Inode/")
count<- read.table("ChrABC-Inode.list",header = T,row.names = 1)
count <- count[rowMeans(count) >1,]
data <- read.table("DF-simple.list",header = T,row.names = 1,sep="\t")
data_group <- as.factor(data$type)#分组信息
group_frame <- as.data.frame(table(data_group))#table转化数据框
names(group_frame) <- c("group","num") #数据框重新命名
group <- as.character(group_frame$group)  #获得分组
num <-  as.character(group_frame$num)  #获得每组的数目
#count保存count数据
#group保存样本分子信息

Batch_Deseq_differnece<-function(exprSet,group,num,save_dir="Alldiffenece",save_dir2="NEW_MA"){
  try({
    ##create a folder 
    save_dir<-paste0(save_dir,"/")
    dir.create(save_dir)
    ## creat a group
    group_list= factor(rep(group,num))
    group_list
    colData=data.frame(row.names = colnames(exprSet),
                       group=group_list)
    
    #dat<-data.frame()
    ## use the Deseq2 to have Diffence analyse
    for (i in 1:length(group)){
      name=unique(group)[i]
      print(name)
      colData$group<-relevel(colData$group,ref=name)
      dds=DESeq2::DESeqDataSetFromMatrix(countData = round(exprSet),#整数化#如果是
                                         colData = colData,
                                         design = ~group) 
      dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
      dds <- DESeq2::DESeq(dds)
      for (j in 2:length(DESeq2::resultsNames(dds))){
        
        resname=DESeq2::resultsNames(dds)[j]
        
        res=DESeq2::results(dds, name=resname)
        
        res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
        res_lfc
        #res=res_lfc
        
        summary(res_lfc)
        summary(res)
        
        dir.create(save_dir)
        
        #write.csv(res,paste0(save_dir,resname,".csv"))
        
        save_dir2=paste0(save_dir2,"/")
        dir.create(save_dir2)
        
        
        
        save_dir_MA=paste0(save_dir2,"/",resname)
        dir.create(save_dir_MA)
        write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
        write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
        
        
        #重新画图
        DEG <- res
        dim(DEG)
        #这里用的是DESeq2的DEG，删掉NA。
        DEG <- na.omit(DEG)
        dim(DEG)
        # 使用基础函数plot绘图
        png(paste0(save_dir_MA,"/",resname,"_MAplot.png"),width=600*3,height=3*600,res=72*3) 
        plot(DEG$log2FoldChange,-log2(DEG$padj))
        dev.off()
        # 确定差异表达倍数,abs表示绝对值
        logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
        # 取前两位小数
        logFC_cutoff <- round(logFC_cutoff, 2)
        logFC_cutoff
        #查看一下logFC_cutoff的值。
        
        # 确定上下调表达基因。
        #方法1：按照logFC_cutoff绘图
        #方法2：按照差异基因的log2FoldChange大于1.5或者其他值绘图。
        DEG$change = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                      ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'STABLE'))
        
        DEG$change = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1.5,
                                      ifelse(DEG$log2FoldChange > 1.5 ,'UP','DOWN'),'STABLE'))
        
        
        
        ##确定需要在火山图上展示的字。
        this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                            '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                            '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
        
        #this_tile <- paste0('Cutoff for logFC is ',round(2.0,1),
        # '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
        #'\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))                    
        
        ###绘图的时候必须是dataframe，所以需要转换一下。
        DEG <- data.frame(DEG)
        write.csv(DEG,paste0(save_dir,resname,"DEG.csv"))
        ## 颜色与分组一一对应
        
        g <- ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj),color=change)) + 
          geom_point(shape = 16, size=2) + 
          theme_set(theme_set(theme_bw(base_size=20))) + 
          xlab("log2 fold change") + 
          ylab("-log10 p-value") +
          #ggtitle( this_tile ) + ##这个可以不要
          theme(plot.title = element_text(size=15,hjust = 0.5)) + 
          theme_classic()+
          scale_colour_manual(values = c('blue','black','red'))
        ggsave(paste0(save_dir_MA,"/",resname,"_MArrput.png"),g) 
        
        png(paste0(save_dir_MA,"/",resname,"_MA.png"),width=600*3,height=3*600,res=72*3) 
        plotMA(res, ylim=c(-3,3),main=paste0(resname," MA"))
        dev.off()
        
        
        png(paste0(save_dir_MA,"/",resname,"_MAlfc.png"),width=600*3,height=3*600,res=72*3) 
        xlim <- c(1,1e5); ylim<-c(-3,3)
        plotMA( res_lfc, xlim=xlim, ylim=ylim, main=paste0(resname,"apeglm"))
        dev.off()
        
      }
      
    }
    
    
  })
}
Batch_Deseq_differnece(count,group=group,num,save_dir = "New",save_dir2="NEW_MA")