#水稻KEGG，GO分析
setwd("D:/bio-information/OS/RNA-SEQ/Metal/DEG/KEGG-GO/")
library(AnnotationHub)
library(ELISAtools)
hub <- AnnotationHub()
hu <- query(hub,"oryza sativa")

Rice <- hub[['AH101067']]

display(Rice)
saveDb(Rice,file="Rice.org")
columns(Rice)
loadDB("Rice.orrdb")
head((keys(Rice,keytype = "ONTOLOGYAL")))


#以上代码，下载好了水稻注释包

library(biomaRt)

a <- read.table("ALL_have_DEG.list")

#在ensemble plants上能看到所有已提交的物种信息
ensembl = useMart(biomart = "plants_mart",host = "http://plants.ensembl.org")
#查看ensemble plants都有哪些物种信息，并设置为该物种信息。
dataset <- listDatasets(mart = ensembl)
head(dataset)
ensembl = useMart(biomart = "plants_mart",host = "http://plants.ensembl.org",dataset="osativa_eg_gene")
#查看该dataset上都有哪些属性，方便后面做添加
attributes <- listAttributes(ensembl)



entr_ID <- getBM(attributes =c("ensembl_gene_id",'entrezgene_id'),
                     filters = "ensembl_gene_id",values = a,mart = ensembl)
gene_name <- getBM(attributes =c("ensembl_gene_id",'external_gene_name',"description"),filters = "ensembl_gene_id",values = a,mart = ensembl)

Go_id <- getBM(attributes =c("ensembl_gene_id",'go_id','goslim_goa_description'),
               filters = "ensembl_gene_id",values = a,mart = ensembl)

