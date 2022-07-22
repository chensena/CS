
setwd("E:/XK4-XK5-RNA-SEQ/New/DEG-XK5_vs_WT/")
install.packages("org.Dfarinosus.eg.db", repos=NULL, type="sources")

library(org.Dfarinosus.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(ggridges)
library(enrichplot)
columns(org.Dfarinosus.eg.db)


DD<-"DOWN-DEG.txt"
DEGs<- read.table(DD, header=TRUE, sep = "\t")
gene_list <- DEGs[,1]

################################################
# 从 OrgDB 提取 Pathway 和基因的对应关系
################################################

pathway2gene <- AnnotationDbi::select(org.Dfarinosus.eg.db, 
                                      keys = keys(org.Dfarinosus.eg.db), 
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

################################################
# 导入 Pathway 与名称对应关系
################################################
load("D:/bash-shell/gene-annotation/orgDB/kegg_info.RData")#储存pathway2name

#KEGG pathway 富集
ekp <- enricher(gene_list, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)

ekp_results <- as.data.frame(ekp)

barplot(ekp, showCategory=20,color="pvalue",
        font.size=10)
dotplot(ekp)

emapplot(ekp)



#########################################################################################
#GO 分析
#########################################################################################

ego <- enrichGO(gene = gene_list,                       #差异基因 vector 
                keyType = "GID",                                   #差异基因的 ID 类型，需要是 OrgDb 支持的 
                OrgDb = org.Dfarinosus.eg.db,                               #对应的OrgDb 
                ont = "MF",                                             #GO 分类名称，CC BP MF 
                pvalueCutoff = 1,                                   #Pvalue 阈值 （pvalue=1指输出所有结果，pvalue=0.05指输出符合要求的结果） 
                qvalueCutoff = 1,                                   #qvalue 阈值 pAdjustMethod = "BH", #Pvalue 矫正方法 
                readable = FALSE)                               #TRUE 则展示SYMBOL，FALSE 则展示原来的ID（选false是因为不是所有gene都有symbol的)

ego_results<-as.data.frame(ego)                          ###生成的ego文件转换成data.frame格式即可。

write.table(ego_results, file = "ego_results.txt", quote = F)                    ###让保存的字符串不用“”引起来
pdf(file = "ego_barplot.pdf")                                                                   ##打开一个PDF文件
barplot(ego, showCategory=20, x = "GeneRatio")                                ##把图画到这个PDF文件里
dev.off()                                                                                                 ##关闭PDF

dotplot(ego)               
emapplot(ego)


for (i in 1:length(ego$Description)){
  
  list<-  strsplit(ego$geneID[i],'/')
  write.csv(list[1],paste0("MF/",ego$Description[i],".csv"),row.names = F,quote=F)
}





setwd("E:/XK4-XK5-RNA-SEQ/New/XK4_vs_WT/")
###GSEA分析
DD<-"DEG.txt"
DEGs<- read.table(DD, header=TRUE, sep = "\t")#读入DESeq2差异表达结果
new_DEGs  <- DEGs[order(-DEGs$log2FoldChange),]#按着log2FoldChange的值排序
gene_list <- new_DEGs$log2FoldChange #提取logfold2的值
names(gene_list) <- new_DEGs$GID#给值命名，构建好genelist
gsego_BP <- gseGO(
                gene_list ,
                ont = "BP",
                OrgDb = org.Dfarinosus.eg.db,
                keyType = "GID",
                exponent = 1,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = FALSE
)
gsego_CC <- gseGO(
  gene_list ,
  ont = "CC",
  OrgDb = org.Dfarinosus.eg.db,
  keyType = "GID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

gsego_MF <- gseGO(
  gene_list ,
  ont = "MF",
  OrgDb = org.Dfarinosus.eg.db,
  keyType = "GID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

#感兴趣的基因集信息提取，这里拆分过氧化相关基因

for (i in 1:length(gsego_CC$core_enrichment)){
    list<-  strsplit(gsego_CC$core_enrichment[i],'/')
    write.csv(list[1],paste0("CC/",gsego_CC$Description[i],".csv"),row.names = F,quote=F)
}

for (i in 1:length(gsego_MF$core_enrichment)){
  list<-  strsplit(gsego_MF$core_enrichment[i],'/')
  write.csv(list[1],paste0("MF/",gsego_MF$Description[i],".csv"),row.names = F,quote=F)
}

for (i in 1:length(gsego_BP$core_enrichment)){
  list<-  strsplit(gsego_BP$core_enrichment[i],'/')
  write.csv(list[1],paste0("BP/",gsego_BP$Description[i],".csv"),row.names = F,quote=F)
}

#GseGO结果可视化

dotplot(gsego_BP,split=".sign")+facet_wrap(~.sign,scales = "free") #换个显示方式,激活，抑制分开显示
dotplot(gsego_CC,split=".sign")+facet_wrap(~.sign,scales = "free") #换个显示方式,激活，抑制分开显示
dotplot(gsego_MF,split=".sign")+facet_wrap(~.sign,scales = "free") #换个显示方式,激活，抑制分开显示


gseaplot2(gsego_BP,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
gseaplot2(gsego_CC,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
gseaplot2(gsego_MF,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p


gseaplot2(gsego_BP,1:10,color="red") # 按第一个做二维码图，并显示p值
gseaplot2(gsego_CC,1:10,color="red") # 按第一个做二维码图，并显示p值
gseaplot2(gsego_MF,1:10,color="red") # 按第一个做二维码图，并显示p


ridgeplot(gsego_BP)#山脊图
ridgeplot(gsego_CC)
ridgeplot(gsego_MF)





gseaplot(gseGO_up,geneSetID = 2,title = gseGO_up$Description[2] )


setwd("D:/bash-shell/gene-annotation/orgDB/")

write.csv(ko2pathway,"ko2pathway.csv")

write.csv(pathway2name,"pathway2name.csv")















#########水稻KEGG-GO基因富集
setwd("D:/bio-information/OS/RNA-SEQ/Metal/DEG/KEGG-GO/KOBAS/")
pathway<-read.table("output_identify_20220721104951.txt",header = T,sep = "\t",stringsAsFactors = F)

pathway$RichFactor<-pathway$Input.number/pathway$Background.number 
kegg <- pathway
library(ggplot2)

p<-ggplot(kegg,aes(x=RichFactor,y=Pathway))


p+geom_point(aes(size=Input.number,color=Corrected.P.Value))+
  scale_color_gradient(low="red",high="green")+
  labs(title="Statistics of Pathway Enrichment",x="Rich factor",y="",color="qvalue",size="Gene_number")+
  theme_bw()

#气泡图

#下面画柱状图

kegg<-read.table("KEGG.list",header = T,sep = "\t",stringsAsFactors = F)

kegg <- kegg[order(kegg$padj),]#根据P值排序

p <- ggplot(data=kegg,aes(x=Description,y=count,fill=padj)) +
  geom_bar(stat="identity") + 
  coord_flip()+
  theme(panel.background=element_rect(fill='transparent',color='blue'),
        axis.text.y=element_text(color="black",size=12))




