library(Biobase)
library(edgeR)
library(DESeq2)
library(GenomicFeatures)
############################先进行筛选 ≥6 reads in at least 20% of samples##########
load("~/RNAseq/alignOut/99sample/99first.RData")
counts.suc <- apply(countsall >= 6,1,sum)
counts.suclist <- counts.suc >=2
counts.res <- countsall[counts.suclist,]
dim(counts.res)
##[1] 18981    99
tpm.res <- tpmall[counts.suclist,]
save(counts.res, tpm.res, phe,file="~/RNAseq/alignOut/99sample/99filter_6reads.RData")
################################导入数据#################################
read.table("~/RNAseq/alignOut/99sample/GeneLengthMeans.txt",header=F,row.name=1) -> genlen
load("~/RNAseq/alignOut/99sample/99filter_6reads.RData")
#提取所需的基因的长度
genlen <- subset(genlen, rownames(genlen) %in% rownames(countsall))
colnames(genlen)<-"EffectiveLength"
#把整理好的用于求表达量的矩阵都保存起来
save(genlen,countsall,tpmall,phe,file="~/RNAseq/alignOut/99sample/99genlen_counts_phe_tpm.RData")
####利用edgeR进行TMM矫正######
len<-t(as.numeric(genlen$EffectiveLength))
all.y<-DGEList(counts=countsall)
all.y1<-calcNormFactors(all.y)
all.y2<-estimateDisp(all.y1)
all.y2<-estimateCommonDisp(all.y2)
all.y2<-estimateTagwiseDisp(all.y2)
norm_counts.table <- t(t(all.y2$pseudo.counts)*(all.y2$samples$norm.factors))
save(norm_counts.table,file="~/RNAseq/alignOut/99sample/NormalizedReadsCountTMM.RData")
############转换为TPM##########
all.y2$genes$Length<-c(len)#
rpkm(all.y2)->all.y3
save(all.y3,file="~/RNAseq/alignOut/99sample/allRPKM_afterTMMRNASeqAll.RData")
apply(all.y3,2,sum)->all.total
alltpm=t((t(all.y3)/all.total)*(10^6))
save(alltpm,file="~/RNAseq/alignOut/99sample/Alltpm.AfterTMMRNASeqAll.RData")
###对低表达量基因进行过滤###################################
tpm.suc <- apply(alltpm>=0.1,1,sum)
tpm.suclist <- tpm.suc>=2
tpm.res <- alltpm[tpm.suclist,]
all.y3[tpm.suclist,]->rpkm.res
save(rpkm.res,file="~/RNAseq/alignOut/99sample/allrpkm.filtered.RNASeqAll.RData")
save(tpm.res,file="~/RNAseq/alignOut/99sample/alltpm.filtered.RNASeqAll.RData")
####把这里的数据送去做deconvolution######
#####做INT转化#######
load("alltpm.filtered.RNASeqAll.RData")
tpm.res.int <- matrix(,nrow(tpm.res),ncol(tpm.res))
for (i in 1:nrow(tpm.res)){
  tpm.res.int[i,] <- qqnorm(tpm.res[i,],plot.it=F)$x
}
tpm.res.int <- data.frame(tpm.res.int)
rownames(tpm.res.int)<-rownames(tpm.res)
colnames(tpm.res.int)<-colnames(tpm.res)
tpm.res.int->tpm.TPMFinal
save(tpm.TPMFinal,file="~/RNAseq/alignOut/99sample/allTPMAfterINT.RData")
#cd /home/shengxin/PEER/AllRNAseq
load("alltpm.filtered.RNASeqAll.RData")
rownames(tpm.res)->humanSymbols
write.table(tpm.res,file="~/RNAseq/alignOut/99sample/BulkCounts.txt",sep="\t",quote=F)