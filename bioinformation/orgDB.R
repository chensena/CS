library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
#顺手设置一下options
options(stringsAsFactors = F)

setwd("D:/bash-shell/gene-annotation/orgDB/")
emapper <- read.table("Dfa.annotations.txt", header=TRUE, sep = "\t",quote = "")
#将空值替换为NA，方便后续使用na.omit()函数提出没有注释到的行
emapper[emapper==""]<-NA

gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()

gene2go = data.frame(GID = character(),
                    GO = character(),
                    EVIDENCE = character())
gos_list <- function(x){
    the_gos <- str_split(x[2],",",simplify=FALSE)[[1]]
    df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                          GO = the_gos,
                          EVIDENCE = rep("IEA", length(the_gos)))
    return(df_temp)
}
gene2gol <- apply(as.matrix(gos),1,gos_list)
gene2gol_df <- do.call(rbind.data.frame, gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)


gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2ko<-na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko),1,gos_list)
gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)

update_kegg <- function(json = "ko00001.json"){
  pathway2name <- tibble(Pathway=character(),Name=character())
  ko2pathway <- tibble(Ko=character(),Pathway=character())
  kegg <- fromJSON(json)
  for(a in seq_along(kegg[["children"]][["children"]])){
    A<-kegg[["children"]][["name"]][[a]]
    for(b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])){
      B<- kegg[["children"]][["name"]][[a]]
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5}","")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["nama"]]
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}}
  save(pathway2name, ko2pathway, file = "kegg_info.RData")}
# 调用函数后在本地创建kegg_info.RData文件，以后只需要载入 "kegg_info.RData"即可
update_kegg()
# 载入kegg_info.RData文件
load(file = "kegg_info.RData")             

gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID,Pathway)%>% na.omit()
write.csv(gene2pathway,"gene2pathway.csv")
write.csv(ko2pathway,"ko2pathway.csv")


tar_id = "862986"
genus = " Dendrocalamus"
species = "farinosus"


gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]
gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
gene_info <- gene_info[!duplicated(gene_info),]
#去重
gene2go <- gene2go[,c(1,2,3)]
str(gene2go)

gene2go$GO <- gsub('["]','',gene2go$GO)
str(gene2go)
gene2ko$Ko <- gsub('["]','',gene2ko$Ko)


#建库
makeOrgPackage(gene_info=gene_info,
              go=gene2go,
              ko=gene2ko,
              pathway=gene2pathway,
              version="1.34.1", #版本，使用？makeOrgPackage，拉到最下面查看
              maintainer = "864306157@qq.com", #修改为你的名字和邮箱
              author = "864306157@qq.com", #修改为你的名字和邮箱
              outputDir = ".", #输出文件位置
              tax_id=tar_id, #你在NCBI上查并定义的id
              genus=genus,
              species=species,
              goTable="go")
              





library(clusterProfiler)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)
#hub <- AnnotationHub::AnnotationHub()
#query(hub,"Theobroma cacao")

#STEP1：自己构建的话，首先
setwd("D:/bash-shell/gene-annotation/orgDB/")

egg_f <- "Dfa.emapper.annotations"
egg <- read.table(egg_f, header=TRUE, sep = "\t",quote="",fill=T)

#egg_f <- "test111"
#egg <- read.table(egg_f, header=TRUE, sep = "\t")
#gene_info_new<-read.table("id.txt" , header=TRUE, sep = "\t")


egg[egg=="-"]<-NA
gene_info <- egg %>% 
  dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()

#STEP3-1：挑出query_name与GO注释信息
gterms <- egg %>%
  dplyr::select(query, GOs) %>% na.omit()

#STEP3-2：我们想得到query_name与GO号的对应信息
# 先构建一个空的数据框(弄好大体的架构，表示其中要有GID =》query_name，GO =》GO号， EVIDENCE =》默认IDA)
# 关于IEA：就是一个标准，除了这个标准以外还有许多。IEA就是表示我们的注释是自动注释，无需人工检查http://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)
# 两种情况下需要用IEA：1. manually constructed mappings between external classification systems and GO terms； 2.automatic transfer of annotation to orthologous gene products.
gene2go <- data.frame(GID = character(),
                      GO = character(),
                      EVIDENCE = character())

# 然后向其中填充：注意到有的query_name对应多个GO，因此我们以GO号为标准，每一行只能有一个GO号，但query_name和Evidence可以重复
for (row in 1:nrow(gterms)) {
  gene_terms <- str_split(gterms[row,"GOs"], ",", simplify = FALSE)[[1]]  
  gene_id <- gterms[row, "query"][[1]]
  tmp <- data_frame(GID = rep(gene_id, length(gene_terms)),
                    GO = gene_terms,
                    EVIDENCE = rep("IEA", length(gene_terms)))
  gene2go <- rbind(gene2go, tmp)
} 

gene2go$GO <- gsub('["]','',gene2go$GO)

#STEP4-1: 挑出query_name与KEGG注释信息
gene2ko <- egg %>%
  dplyr::select(GID = query, KO = KEGG_ko) %>%
  na.omit()

#STEP4-2: 得到pathway2name, ko2pathway
# 需要下载 json文件(这是是经常更新的)
# https://www.genome.jp/kegg-bin/get_htext?ko00001
# 代码来自：http://www.genek.tv/course/225/task/4861/show

if(F){
  # 需要下载 json文件(这是是经常更新的)
  # https://www.genome.jp/kegg-bin/get_htext?ko00001
  # 代码来自：http://www.genek.tv/course/225/task/4861/show
  library(jsonlite)
  library(purrr)
  library(RCurl)
  
  update_kegg <- function(json = "ko00001.json") {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    
    kegg <- fromJSON(json)
    
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
    save(pathway2name, ko2pathway, gene2go, gene2ko, gene_info, file = "kegg_info.RData")
  }
  

  update_kegg(json = "ko00001.json")
  
  save(pathway2name, ko2pathway, gene2go, gene2ko, gene_info, file = "kegg_info.RData")

  
}

#STEP5: 利用GO将gene与pathway联系起来，然后挑出query_name与pathway注释信息
load(file = "kegg_info.RData")

ko2pathway$Ko <- as.character(ko2pathway$Ko)
ko2pathway$Pathway <- as.character(ko2pathway$Pathway)

gene2ko$GID <- as.character(gene2ko$GID)
gene2ko$KO <- as.character(gene2ko$KO)
str(gene2ko)
str(ko2pathway)

gene2go[gene2go=="-"] <- NA
gene2ko[gene2ko=="-"] <- NA

gene2ko$Ko <- gsub('["]','',gene2ko$Ko)
gene2go$GO <- gsub('["]','',gene2go$GO)
colnames(gene2ko) <- c("GID","Ko")

colnames(ko2pathway) <- c("Ko","Pathway")
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub('[Ko:]','',gene2ko$Ko)

#保存以上的分析结果
write.csv(gene2go,"gene2go.csv")
write.csv(gene2ko,"gene2ko.csv")
write.csv(gene_info,"gene_info.csv")
write.csv(ko2pathway,"ko2pathway.csv")

gene2go<-read.csv(file="gene2go.csv",header=T)
gene2ko<-read.csv(file="gene2ko.csv",header=T)
gene_info<-read.csv(file="gene_info.csv",header=T)
ko2pathway<-read.csv(file="ko2pathway.csv",header=T)



gene2pathway <- merge(gene2ko,ko2pathway,by="Ko")

#gene2pathway <- left_join(gene2ko,ko2pathway,by= "Ko") %>% dplyr::select(GID,Pathway) %>% na.omit()
write.csv(gene2pathway,"gene2pathway.csv")

gene2ko <- as.data.frame(gene2ko)

write.csv(gene2pathway,"gene2pathway.csv")
library(AnnotationForge)  


#STEP6： 制作自己的Orgdb
# 查询物种的Taxonomy，例如要查sesame

#gene2go <- unique(gene2go)

tar_id = "862986"
genus = " Dendrocalamus"
species = "farinosus"

#gene2go<-gene2go[!duplicated(gene2go),]
#gene2ko<-gene2ko[!duplicated(gene2ko),]
#gene2pathway<-gene2pathway[!duplicated(gene2pathway),]

makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               maintainer = "864306157@qq.com",
               author = "SWUST,Bamboo Research Institute,chensen",
               pathway=gene2pathway,
               version="0.0.1",
               outputDir = ".",
               tax_id=tar_id,
               genus=genus,
               species=species,
               goTable="go")

ricenew.orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")




#抄袭第三个代码
library(tidyr)
setwd("D:/bash-shell/gene-annotation/orgDB/")
##名称由makeOrgPackage函数规定请不要更改
emapper <- read.delim("Dfa.emapper.annotations")  
emapper[emapper=="-"] <- NA
emapper <- dplyr::select(emapper,GID=query,Gene_Symbol=Preferred_name, 
                GO=GOs,KO=KEGG_ko,Pathway =KEGG_Pathway, 
                OG=eggNOG_OGs,Gene_Name =seed_ortholog)
#GID=query,这里的query就是csv文件中的表头信息
#这里共提取了gene_info,  gene2go,gene2ko,gene2pathway,gene2symbol
#你用哪些信息就可以进行相应的增减。
#如果你只是做go富集分析，其实gene_info和gene2go就足够了
#gene2symbol在这里感觉纯粹是凑数，这是参考引文2弄的。
#gene2ko，gene2pathway是用来做kegg富集分析的。

#提取GID与Gene_Name信息，参考1是将X.4作为Gene_Name，应该是自己定义的信息，这里将seed_ortholog作为Gene_Name信息，你可以根据实际再调整。
gene_info <- dplyr::select(emapper,GID,Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))

#提取GID与GO信息，组成goTable，建库时的goTable需要三列（GID, GO和EVIDENCE），少一列，就会出错.
gene2go <- dplyr::select(emapper,GID,GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  dplyr::filter(GO!="-",!is.na(GO))%>%   #这是只提取有GO注释信息的行，判断的标准时GO信息不是NA，这也就是为什么前面必须将“-”替换为NA，不替换就无法进行有效过滤。
  mutate(EVIDENCE = 'A')     #硬生生加了1列EVIDENCE，全部赋值A,凑数的。
dim(gene2go)    #查看数据维度。
#[1] 1523399       3

#提取GID与KO信息，这里只有2列信息
gene2ko<- dplyr::select(emapper,GID,KO) %>%
  separate_rows(KO, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(KO))
dim(gene2ko)
#[1] 30530     2

#提取GID与Pathway信息，这里只有2列信息
gene2pathway<- dplyr::select(emapper,GID,Pathway) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(Pathway))
dim(gene2pathway)
#[1] 143056      2

#提取GID与Gene_Symbol信息，Gene_Symbol是Preferred_name信息，这里只有2列信息
gene2symbol<- dplyr::select(emapper,GID,Gene_Symbol) %>%
  dplyr::filter(!is.na(Gene_Symbol))
dim(gene2symbol)
#[1] 4561    2


tar_id = "862986"
genus = "Dendrocalamus"
species = "farinosus"



library(AnnotationForge)  #如果没有就安装一下。
AnnotationForge::makeOrgPackage(gene_info=gene_info,
                                go=gene2go,
                                ko=gene2ko,
                                pathway=gene2pathway,
                                symbol=gene2symbol,
                                maintainer='CS<864306157@qq.com>',
                                author='ChenSen,SWUST',
                                version="0.1" ,
                                outputDir=".", 
                                tax_id=tar_id,
                                genus=genus,
                                species=species,
                                goTable="go")




