library(multcomp)


#多重方差检验
#只需要一个输入文件
setwd("E:/Sub_chr_E/OS-DFA-DFB-DFC/")
p <- read.table('DFA-DFB-DFC.list',header = TRUE,sep = "\t")

p$Group <- as.factor(p$Group)
fit<-aov(Ka.Ks~Group,data=p)
TukeyHSD(fit)

par(mar=c(5,4,6,4))

tuk<-glht(fit,linfct=mcp(Group="Tukey"))
plot(cld(tuk,level=0.05),col="lightgreen")



#T检验
library("ggpubr")
setwd("E:/Sub_chr_E/expersion/ABC-expersion/Inode/New/xianzhu/")
data = read.table("B_vs_C.list",header=T)
com1 <- compare_means( lg2Fc~ group , data = data, method = "t.test")


p <- ggboxplot(data, x="group", y = "lg2Fc", color = "group", palette = "jco", add = "jitter",  short.panel.labs = FALSE) +

# 添加p值
stat_compare_means(method = "t.test",label.y=100) +
# 显示p值但不显示方法
stat_compare_means(aes(label = ..p.format..),method = "t.test",label.x = 1.5)+

# 只显示显著性水平
stat_compare_means(aes(label = ..p.signif..),method = "t.test",label.x = 1.5)

