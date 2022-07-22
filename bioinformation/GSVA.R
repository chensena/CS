library(GSVA)
library(org.Dfarinosus.eg.db) #导入梁山慈竹的注释包
columns(org.Dfarinosus.eg.db)
library('GSEABase')


#制作梁山慈竹数据集
rm(list=ls()) #清除环境变量
setwd("D:/bio-information/DF/GSVA/")
#GO数据集
goannot <- AnnotationDbi::select(org.Dfarinosus.eg.db, keys=keys(org.Dfarinosus.eg.db), columns=c("GO","EVIDENCE","ONTOLOGY"))
genesbygo <- split(goannot$GID, goannot$GO)
#KO数据集合
load("D:/bash-shell/gene-annotation/orgDB/kegg_info.RData")#储存pathway2name
koannot<- AnnotationDbi::select(org.Dfarinosus.eg.db, keys=keys(org.Dfarinosus.eg.db), columns=c("KO","Pathway"))
genesbyko <- split(koannot$GID, koannot$KO)

#读入表达矩阵
datExpr <- read.table("All.DEG_final.txt",header=T,row.names = 1)
#表达矩阵取log值，用于做归一化
datExpr <- log2(datExpr+1)
datExpr <- as.matrix(datExpr)

#富集到KEGG ko
gsva.es <- gsva(datExpr,genesbyko,method='gsva',kcdf='Gaussian',verbose=T, 
                parallel.sz = parallel::detectCores())#gsva计算







gsva.df <- as.data.frame(gsva.es)
gsva2 <- cbind(rownames(gsva.df),gsva.df)#将KO加入数据框中


#富集热图
pheatmap::pheatmap(gsva.es)

#ko号转换为name
gsva2 <- cbind(rownames(gsva.df),gsva.df)#将KO加入数据框中
colnames(gsva2)[1] <- "Ko"  #给新加的列取名
gsva2 <- as.data.frame(gsva2)
gsva2$KO <- gsub('[ko:]','',gsva2$KO)#去除KO列的“Ko:”
#ko2pathWay  将pathway 添加到ko后
gsva3 <- merge(gsva2,ko2pathway,by="Ko")
#pathway2name  #将pathway name 添加到pathway后面
gsva4 <- merge(gsva3,pathway2name,by="Pathway")


library(limma)
#limmad多重差异表达分析
#以下为设置分组信息
data <- read.table("DF-sample.list",header = T,row.names = 1,sep="\t")
data_group <- as.factor(data$type)#分组信息
group_frame <- as.data.frame(table(data_group))#table转化数据框
names(group_frame) <- c("group","num") #数据框重新命名
group <- as.character(group_frame$group)  #获得分组
num <-  as.character(group_frame$num)  #获得每组的数目
condition <- factor(rep(group,num),levels=group)
table(condition)#多此一举，实际上第55行就可以了

#设定差异比较矩阵 **这里注意了，经常绕不清楚的地方来了**
# 这是需要声明差异比较矩阵的方法
design <- model.matrix(~0+factor(condition))
colnames(design) <- levels(factor(condition))
# Tunor VS Normal
compare <- makeContrasts(Root - MLeaf, levels=design)
fit <- lmFit(gsva.es, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Dif <- topTable(fit3, coef=1, number=200)
head(Diff)


#对GSVA差异分析结果进行热图可视化

padj_cutoff=0.05
log2FC_cutoff=log2(2)
degs <- Dif
keep <- rownames(degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2FC_cutoff, ])
length(keep)
dat <- gsva.es[keep[1:50],] #选取前50进行展示

pheatmap::pheatmap(dat, 
                   fontsize_row = 8,
                   height = 10,
                   width=18,
                   annotation_col = gl,
                   show_colnames = F,
                   show_rownames = T,
                   filename = paste0('GSVA_go_heatmap.pdf'))

#发散条形图绘制

dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
install.packages("ggprism")
install.packages("ggthemes")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签
# 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标


#火山图

degs$significance  <- as.factor(ifelse(degs$adj.P.Val < padj_cutoff & abs(degs$logFC) > log2FC_cutoff,
                                       ifelse(degs$logFC > log2FC_cutoff ,'UP','DOWN'),'NOT'))

this_title <- paste0(' Up :  ',nrow(degs[degs$significance =='UP',]) ,
                     '\n Down : ',nrow(degs[degs$significance =='DOWN',]),
                     '\n adj.P.Val <= ',padj_cutoff,
                     '\n FoldChange >= ',round(2^log2FC_cutoff,3))

g <- ggplot(data=degs, 
            aes(x=logFC, y=-log10(adj.P.Val),
                color=significance)) +
  #点和背景
  geom_point(alpha=0.4, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( adj.P.Val )") +
  #标题文本
  ggtitle( this_title ) +
  #分区颜色                   
  scale_colour_manual(values = c('blue','grey','red'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(),
        legend.position="right")

ggsave(g,filename = 'GSVA_go_volcano_padj.pdf',width =8,height =7.5)

