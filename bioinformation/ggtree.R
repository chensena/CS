library(ggtree)
library(treeio)
library(tidyverse)
library(ggsci)

setwd("D:/bio-information/DF/genefamily/laccase/PHY/OS-ZM-Dfa/best-Model/")
tree <- read.iqtree("D:\\bio-information\\SUBCHrsome\\DF-BAM-OS-OLA-GAN-BAM\\Single_Copy\\new\\all.muscle.trim.fa.treefile")

tree <- read.tree("D:/bio-information/DF/genefamily/laccase/PHY/OS-ZM-Dfa/best-Model/Os-Zm-Dfa-LAC.muscle.trim.sim.fas.contree")
#读取进化树文件

group <- read.table("group.list",header = T,row.names=1)
groupinfo <- split(row.names(group),group$group)
#ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
tree <- groupOTU(tree,groupinfo)

p <- ggtree(tree,branch.length = "none",linetype = 1,size=0.8,ladderize = T,aes(color=group))



p <- groupOTU(tree,groupinfo) %>%
  ggtree(branch.length = "none",
         linetype = 1,layout="circular",size=0.8,ladderize = T,aes(color=group))
P1 <- p + #geom_text(aes(label=node),hjust=-1)+
  scale_color_nejm() + labs(color="") +
  geom_tiplab(aes(color=group),size = 5,hjust = 0) +
  xlim(NA,22) +
  geom_text2(aes(subset=!isTip,label=node),hjust=-.3,color="red")
  
  geom_highlight(node=93,fill = "#FFB6C1",type="rect",alpha=0.6)+
  geom_highlight(node=87,fill = "#E6E6FA",type="rect",alpha=0.6)+
  geom_highlight(node=69,fill = "#E0FFFF",type="rect",alpha=0.6) +
  geom_highlight(node=2,fill = "#E0FFFF",type="rect",alpha=0.6) +
  geom_cladelabel(node=17, label = "chrC", align = T, offset = 3.3, color = "red", barsize = 2,fontsize = 6) +
  geom_cladelabel(node=20, label = "chrA", align = T, offset = 3.3, color = "pink", barsize = 2,fontsize = 6) +
  geom_cladelabel(node=22, label = "chrB", align = T, offset = 3.3, color = "blue",barsize =2,fontsize = 6) +
  geom_cladelabel(node=2, label = "H", align = T, offset = 3.3, color = "green",barsize = 2,fontsize = 6)

P2 <- ggtree::rotate(P1,60) 
P3 <- ggtree::rotate(P2,59)+
  geom_highlight(node=93,fill = "#FFB6C1",type="rect",alpha=0.6)+
  geom_highlight(node=87,fill = "#E6E6FA",type="rect",alpha=0.6)+
  geom_highlight(node=69,fill = "#E0FFFF",type="rect",alpha=0.6) +
  geom_highlight(node=79,fill = "#E0FFFF",type="rect",alpha=0.6)+
  geom_highlight(node=61,fill = "#AFEEEE",type="rect",alpha=0.6)+
  geom_nodelab(hjust=-0.05,size=4)
  


