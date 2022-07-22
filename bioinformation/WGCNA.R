library(WGCNA)
library(reshape2)
library(stringr)
#Z载入WGCNA程序包
###
setwd("E:/XK4-XK5-RNA-SEQ/")
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads(nThreads=10)
exprMat <- "All.DEG_final.xls"
type = "signed"

dataExpr <- read.table(exprMat, sep='\t', header=T, row.names = 1,
                       quote="", comment="", check.names=F)
#读入差异表达矩阵

#数据筛选

m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25),)[2],0.01)),]
###筛选出52445个基因

####下面依旧是一种筛选方法



################################

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))


## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}



nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#经过检测，暂时没有发现离群样本
#保存完毕，文件保存位置 “D:\\bio-information\\WGCNA”

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)



  
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")


# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "signed", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "signed", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "signed", 7, 14),
                               ifelse(type == "signed", 6, 12))      
                 )
  )
}

#power = 8

##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
power=sft$powerEstimate
power=18#洋黄竹转录组表达矩阵，没有筛选到合适的软阈值，用经验power值
power <- as.numeric(power)
#power = 12#未找到
#robustY = ifelse(corType=="pearson",T,F)
#corType可以不用设置


#如果net报错
dataExpr[] <- lapply(dataExpr, as.numeric)

#power 软阈值 如果没有找到合适的软阈值，unsigned 选6，singned选12
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, 
                        loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)


#she'zhi
##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will use 47 parallel threads.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file WGCNA/LiverFemaleClean.txt.tom-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 3 genes from module 1 because their KME is too low.
##      ..removing 5 genes from module 12 because their KME is too low.
##      ..removing 1 genes from module 14 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.25
##        Calculating new MEs...

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)
##绘画结果展示### open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#保存模型结果


#0     1     2     3     4     5     6     7     8     9    10    11    12    13 
#1286 11987 10552  5758  4359  3959  3733  2528  2026  1437  1359   591   501   488 
#14    15    16    17    18    19    20    21    22    23    24 
#351   309   297   248   172   149   107    95    57    53    43 


##结果保存###
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
x <- table(moduleColors)
write.csv(x,"模块基因数目.csv")
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,net,sft,TOM,
     file = "All-color-02-networkConstruction-auto.RData")

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
####导入上一步分析的数据
l <- load("All-color-02-networkConstruction-auto.RData")
k <- load("All.DEG_final.xls.tom-block.1.RData")
j <- load("性状模块关联-01-dataInput.RData")


##########################


MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)




###导出网络到Cytoscape#### Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(dataExpr, power = power)
# Read in the annotation file# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
#modules = c("blue","brown","green","magenta","pink","red","turquoise","yellow")
# Select module probes选择模块探测




#表型数据拟采用木质素，纤维素含量，还有各种激素的含量，时间戳，统统加上
#2022.5.31 模块与性状联系"
trait <- "Sample.list"
traitData <- read.table(trait, sep='\t', header=T, row.names=1,
                                   check.names=FALSE, comment='',quote="")



corType = "pearson"#严重依赖二元变量
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)


if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
               cex.lab = 1,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.7, zlim = c(-1,1),
               main = paste("Module-trait relationships"))


##表型相关性结束

####2022,6,4,性状与表型关联研究

nGenes = ncol(dataExpr);
nSamples = nrow(dataExpr);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData, use = "p",method="spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# 用热图的形式展示相关系数

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#该分析确定了几个重要的模块-特征关联。我们将体重作为感兴趣的特征来研究。

















#以下代码作废
moduleColors = mergedColors 
table(moduleColors)
colorOrder= c( "grey", standardColors( 50)) 
moduleLabels= match(moduleColors, colorOrder)- 1 
MEs= mergedMEs

nGenes= ncol(dataExpr) 
nSamples= nrow(dataExpr) 
moduleTraitCor= cor(MEs, datTraits, use = "p",method="spearman") 
moduleTraitPvalue= corPvalueStudent(moduleTraitCor, nSamples)


textMatrix= paste(signif(moduleTraitCor, 2), "n(", signif(moduleTraitPvalue,1), ")", sep = "") 
dim(textMatrix)= dim(moduleTraitCor) 
par(mar= c(5, 10, 3, 3)) 
labeledHeatmap(Matrix= moduleTraitCor, 
               xLabels= names(datTraits), 
               yLabels= names(MEs),
               ySymbols= names(MEs), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE,
               cex.text= 0.6, zlim= c(-1,1), 
               main= paste("Module-trait relationships"))

####################性状与表型关联





#####################################################共表达基因统计
traitData <- read.table("D:/bio-information/WGCNA/ALL-SAMPLE-WGCNA/sample.list",header = T)
dim(traitData)  #每行是一个样本，每列是一种信息
names(traitData)
#删除我们不需要的数据

#####本脚本部进行下面四行的删除操作,非数值型性可以根据以下操作剔除
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];  #只保留数值型数据
dim(allTraits)
names(allTraits)
########


Samples = rownames(dataExpr,);
traitRows = match(Samples, traitData$ID);
datTraits = traitData [traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage()


sampleTree2 = hclust(dist(dataExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


save(dataExpr, datTraits, file = "性状模块关联-01-dataInput.RData")





module = "organge"

geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
modNames = substring(names(MEs), 3)
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2@@@@@@














y <- as.data.frame(table(moduleColors))#table转化数据框
names(y) <- c("color","num") #数据框重新命名
color <- as.character(y$color)  #获得分组

for (i in color){

modules = i
probes = names(dataExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,                               
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
}


## 可视化基因网络## 
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = power)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap

plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot

diag(plotTOM) = NA
# Call the plot function#sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#随便选取1000个基因来可视化
nSelect = 500
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There''s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
sizeGrWindow(7, 6) 
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)


#模块合并
MEDissThres = 0.4
merge_modules = mergeCloseModules(five, moduleColors, cutHeight = MEDissThres, verbose
                                  = 3)

