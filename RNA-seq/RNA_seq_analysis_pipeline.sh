######0.QC#######
trim_galore --fastqc --retain_unpaired --paired /data/center/11_LN_Chen/5_rnaseq/01.RawData/LR1_1.fq.gz /data/center/11_LN_Chen/5_rnaseq/01.RawData/LR1_2.fq.gz -o /public/home/zhaohuanhuan/RNAseq/test/chuli
###1.Index######
cd ~/lunix_lesson/rnaseq/raw_data
mkdir -p /n/scratch2/username/chr1_hg38_index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/hg19/ --genomeFastaFiles /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/GRCh37.p13.genome.fa --sjdbGTFfile /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/gencode.v19.annotation.gtf --sjdbOverhang 149 ####sjdbOverhang测序长度最大值减1
###2.Align####
mkdir ../results/STAR #
STAR --runMode alignReads --runThreadN 16 \  
--genomeDir /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/hg19 \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--readFilesIn /public/home/zhaohuanhuan/RNAseq/test/chuli/LR1_1_val_1.fq.gz /public/home/zhaohuanhuan/RNAseq/test/chuli/LR1_2_val_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix LN1 \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outSAMunmapped Within KeepPairs \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType WithinBAM SoftClip \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outSAMattrRGline ID:rg1 SM:sm1
##3.MarkDuplicates##
cat /public/home/zhaohuanhuan/RNAseq/alignOut/RSEM/samplelist.txt | while read i 
do
python3 -u /public/home/zhaohuanhuan/eQTL/gtex-pipeline-master/rnaseq/src/run_MarkDuplicates.py /public/home/zhaohuanhuan/RNAseq/alignOut/STAR/${i}Aligned.sortedByCoord.out.bam ${i}
done
#####4.RSEM###########
mkdir rsem
cd rsem
#####4.1RSEM-Index###
rsem-prepare-reference --gtf /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/gencode.v19.annotation.gtf \  ####gtf文件
--star \  ###STAR
-p 16 \
/public/home/zhaohuanhuan/RNAseq/reference/index/STAR/GRCh37.p13.genome.fa \  
/public/home/zhaohuanhuan/RNAseq/reference/RSEM/hg19 
####4.2transcript quantification###
cd /public/home/zhaohuanhuan/RNAseq/test/RSEM
rsem-calculate-expression --num-threads 16 \
--fragment-length-max 1000 \
--no-bam-output \
--paired-end \
--estimate-rspd \
--forward-prob 0 \
--bam /public/home/zhaohuanhuan/RNAseq/test/chuli/LN1Aligned.toTranscriptome.out.bam \
/public/home/zhaohuanhuan/RNAseq/reference/RSEM/hg19 \ 
LN1  ####
##4.3###conduct counts-matrix using RSEM results
rsem-generate-data-matrix *.genes.results > /public/home/zhaohuanhuan/RNAseq/alignOut/deseq_out/output.matrix
###
a=`ls *.genes.results | tr "\n" " "`
paste $a > /public/home/zhaohuanhuan/RNAseq/alignOut/deseq_out/alloutput.txt
###导入R#####
read.table("~/RNAseq/alignOut/deseq_out/alloutput.txt",header=T,row.name=1)->all
tpmall=all[c(grep("TPM",colnames(all)))] ####提取colnames包含TPM的列
###用正确顺序的countsall的列名替换现有列名##############
read.table("~/RNAseq/alignOut/deseq_out_STAR100/combinesamplelist.txt",header=F,row.name=1) -> sample  ###header=T会少一行
a <- rownames(sample)
a -> colnames(tpmall)
write.table(tpmall,file="~/RNAseq/alignOut/deseq_out/tpmall.txt",quote=F,sep="\t",row.name=T)
read.table("~/RNAseq/alignOut/deseq_out/output.matrix",header=T,row.name=1)->countsall
a -> colnames(countsall)
#########666666666666666666666666666提取effective_length用于计算基因长度####################
read.table("~/RNAseq/alignOut/deseq_out/alloutput.txt",header=T,row.name=1)->all
length=all[c(grep("effective_length",colnames(all)))] ####提取colnames包含effective_length的列
write.table(length,file="~/RNAseq/alignOut/99sample/AllMergedGenesELength.txt",quote=F,sep="\t",row.name=T)
###############取行平均值#########此处未去除最后一个基因
read.table("AllMergedGenesELength.txt",header=T,row.name=1)->genes
last<-dim(genes)[1]
genes[-last,]->genes4
genes4->genes2
a<-matrix(as.numeric(as.character(unlist(genes4))),nrow=dim(genes2)[1])
rownames(a)<-rownames(genes2)
colnames(a)<-colnames(genes2)
glen<-rowMeans(a)
write.table(glen,file="GeneLengthMeans_removelast.txt",col.name=F,quote=F,sep="\t")
#########77777777777777777777777777777############QC
library(DESeq2)
library(proxy)
load("99First.RData")
log2(tpmall+1)->tpm1
prcomp(tpm1,scale=TRUE,center=TRUE)->lgPCA
pdf(file="tpmPCA.pdf")
plot(lgPCA)
dev.off()
countsall2 <- data.frame(apply(countsall,2,function(x)as.numeric(as.character(x))))
round(countsall2) -> countsall3
rownames(countsall3) <- rownames(countsall)
log2(countsall3+1) -> countsall3_1
rowMeans(countsall3_1)->GeneMeans
#colMeans(countsall3)->SampleMeans
density(GeneMeans)->d 
pdf("Distribution of gene mean raw reads logcounts-all.pdf")
#hist(GeneMeans,freq=F,breaks=100)
plot(d)
dev.off()
####
rownames(phe) <- phe$id
dds<-DESeqDataSetFromMatrix(countsall3,DataFrame(phe),design=~id)
dds<-estimateSizeFactors(dds)
normalized_counts<-counts(dds,normalized=TRUE)
rowMeans(normalized_counts)->GeneMeans
density(GeneMeans)->d
pdf("Distribution of gene mean normalized reads.pdf")
plot(d)
dev.off()
log2(1+normalized_counts)->final
rowMeans(final)->GeneMeans
density(GeneMeans)->d
pdf("Distribution of log transformmed gene mean normalized reads counts.pdf")
plot(d)
dev.off()
###########keep genes with a log-transformed value >1 in >2 of samples.
Ssuc<-final>1 
##keepGenes_10<-rowMeans(Ssuc)>0.09615385 
####FALSE  TRUE 
####36012 21808 
keepGenes_2<-rowMeans(Ssuc)>0.02020202 
####FALSE  TRUE 
####31330 26490 
final[keepGenes,]->filtered
rowMeans(filtered)->GeneMeans
density(GeneMeans)->d
pdf("Distribution of log transformmed gene mean normalized reads counts(filtered low exp).pdf")
plot(d)
dev.off()
save(filtered,file="~/RNAseq/alignOut/99sample/filtered99.RData")
scled<-scale(filtered,center=TRUE,scale=TRUE)
save(scled,file="~/RNAseq/alignOut/99sample/scaled99.RData")
rowMeans(scled)->GeneMeans
density(GeneMeans)->d
pdf("Distribution of log transformmed gene mean normalized reads counts(scaled).pdf")
plot(d)
dev.off()
cosine = function(a,b) { 
     len = (sqrt(a %*% a)*sqrt(b %*% b)); 
     if (len == 0) { 
          0; 
     } else { 
          (a %*% b)/len; 
     } 
}
t(scled)->t
hc<-hclust(dist(t,method="cosine"),method="average")
save(hc,file="~/RNAseq/alignOut/deseq_out/cluster104.RData")
pdf("~/RNAseq/alignOut/deseq_out/sampleClusters.pdf",height=10,width=15)
plot(hc)
dev.off()
library("ClassDiscovery")
spca3<-SamplePCA(filtered,usecor=TRUE,center=TRUE)
save(spca3,file="~/RNAseq/alignOut/deseq_out/pcaByCountsNormalized.RData")
head(round(cumsum(spca3@variances)/sum(spca3@variances),digits=3),35)
#[1] 0.164 0.234 0.292 0.338 0.369 0.399 0.425 0.447 0.466 0.482 0.497 0.511
#[13] 0.524 0.536 0.547 0.557 0.567 0.576 0.584 0.593 0.601 0.609 0.617 0.624
#[25] 0.632 0.639 0.646 0.652 0.659 0.666 0.672 0.678 0.685 0.691 0.697
maha15<-mahalanobisQC(spca3,35)
p.adjust(maha15$p.value,"bonferroni")->fdr
data.frame(fdr)->fdr
rownames(maha15)->rownames(fdr)
fdr$color<-"blue"
rownames(fdr[fdr$fdr<0.05,])->outliers
fdr[outliers,]$color<-"red"
prcomp(filtered,scale=TRUE,center=TRUE)->allpca2
pdf("~/RNAseq/alignOut/deseq_out/fdrFailedOutlier.pdf")
plot(allpca2$rotation,col=fdr$color)
dev.off()
length(outliers)
#18
save(allpca2,file="~/RNAseq/alignOut/deseq_out/pcaForPlots.RData")
result2=cutree(hc,k=5) ####
data.frame(result2)->output
as.factor(output$result2)->output$result2
summary(output$result2)
#1   2   3   4   5 
#100   1   1   1   1 
output$color<-"blue"
#rownames(output[output$result2=="2",])->cluster2
rownames(output[output$result2 %in% c("2","3","4","5"),])->cluster2
intersect(cluster2,outliers)->overlap
save(cluster2,file="~/RNAseq/alignOut/deseq_out/Outliers.RData")
write.table(cluster2,file="~/RNAseq/alignOut/deseq_out/Outliersk5.txt",sep="\n",quote=F,row.name=F,col.name=F)
output[cluster2,]$color<-"red"
pdf("~/RNAseq/alignOut/deseq_out/outlierRemovedStep_k5.pdf")
plot(allpca2$rotation,col=output$color)
dev.off()
#########8888888888888888888888筛选基因，生成表达量文件####################
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
save(norm_counts.table,file="NormalizedReadsCountTMM.RData")
############转换为TPM##########
all.y2$genes$Length<-c(len) ######c(len)不能是len
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
write.table(humanSymbols,sep="\n",quote=F,row.name=F,col.name=F,file="~/RNAseq/alignOut/99sample/MyHumanGenes.txt")
###########matrix eQTL需要的基因表达量文件，列是样本，行是基因
write.table(tpm.TPMFinal,file="~/RNAseq/alignOut/99sample/expression/expressionINT.txt",sep="\t",quote=F)