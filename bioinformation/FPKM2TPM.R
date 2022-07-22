#FPKM->TPM
setwd("D:\\bash-shell/rna-seq/")
      
expMatrix<-read.table("D:\\bash-shell/rna-seq/All_gene_fpkm.list",header = T,row.names = 1)


expMatrix<-read.table("D:\\bio-information/MA/rna_seq/SRR/out/SAM/Dla_All_Count.txt",header = T,row.names = NULL)
row.names(mexpMatrix<-make.names(expMatrix[,1],TRUE))
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  fpkm <- as.numeric(fpkm)
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

tpms <- apply(expMatrix,2,fpkmToTpm)
write.table(tpms,'Df_tpm.list',col.names=T,row.names=T,quote=F,sep='\t')
tpms[1:3,]