library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyverse)
library(data.table)
library(sctransform)
#library(clustree)
library(magrittr)
library(RColorBrewer)
library(edgeR)
#library(plyr) -> wilcoxon testing line does not work when this is loaded after dplyr
library(corrplot)
library(pheatmap)
library(ComplexHeatmap)
library(scales)
library(viridis)
library(circlize)
library(ggrepel)
out.path <- "/public/home/zhaohuanhuan/scRNA/B/10HC11LN/DEG/data/"
######22222222对全部细胞的DEG计算##
load("3score_labled15.8immune.combinedCelltypeGroup_celltype_ABC.RData")

meta_data <- immune.combined@meta.data

# 提取 barcode（行名）和 orig.ident 列
barcodes <- rownames(meta_data)
samples <- meta_data$orig.ident

# 将 orig.ident 和 barcode 连接起来生成新的列名
new_colnames <- paste(samples, barcodes, sep = "_")

raw.counts <- immune.combined@assays$RNA@counts
# 将 raw.counts 的列名更新为新生成的列名
colnames(raw.counts) <- new_colnames
raw.counts <- as.data.frame(raw.counts)
raw.counts$hc2 <- rowSums(raw.counts[,grep("hc2", names(raw.counts))])
raw.counts$hc3 <- rowSums(raw.counts[,grep("hc3", names(raw.counts))])
raw.counts$hc4 <- rowSums(raw.counts[,grep("hc4", names(raw.counts))])
raw.counts$hc5 <- rowSums(raw.counts[,grep("hc5", names(raw.counts))])
raw.counts$hc6 <- rowSums(raw.counts[,grep("hc6", names(raw.counts))])
raw.counts$hc7 <- rowSums(raw.counts[,grep("hc7", names(raw.counts))])
raw.counts$hc8 <- rowSums(raw.counts[,grep("hc8", names(raw.counts))])
raw.counts$hc9 <- rowSums(raw.counts[,grep("hc9", names(raw.counts))])
raw.counts$hc10 <- rowSums(raw.counts[,grep("hc10", names(raw.counts))])
raw.counts$hc11 <- rowSums(raw.counts[,grep("hc11", names(raw.counts))])
raw.counts$pbmc1k <- rowSums(raw.counts[,grep("pbmc1k", names(raw.counts))])
raw.counts$pbmc2k <- rowSums(raw.counts[,grep("pbmc2k", names(raw.counts))])
raw.counts$pbmc3k <- rowSums(raw.counts[,grep("pbmc3k", names(raw.counts))])
raw.counts$pbmc4k <- rowSums(raw.counts[,grep("pbmc4k", names(raw.counts))])
raw.counts$pbmc5k <- rowSums(raw.counts[,grep("pbmc5k", names(raw.counts))])
raw.counts$pbmc6k <- rowSums(raw.counts[,grep("pbmc6k", names(raw.counts))])
raw.counts$pbmc7k <- rowSums(raw.counts[,grep("pbmc7k", names(raw.counts))])
raw.counts$pbmc8k <- rowSums(raw.counts[,grep("pbmc8k", names(raw.counts))])
raw.counts$pbmc9k <- rowSums(raw.counts[,grep("pbmc9k", names(raw.counts))])
raw.counts$pbmc10k <- rowSums(raw.counts[,grep("pbmc10k", names(raw.counts))])
raw.counts$pbmc11k <- rowSums(raw.counts[,grep("pbmc11k", names(raw.counts))])
raw.sums <- raw.counts[,237865:237885]
write.csv(raw.sums, file = paste0(out.path, "immune.combined_sample_sum_counts.csv"), row.names = TRUE)
####
cts <- read.csv("~/scRNA/B/10HC11LN/DEG/data/immune.combined_sample_sum_counts.csv", row.names="X")
clinical_data <- read.table("~/scRNA/B/10HC11LN/DEG/data/clinical_data.txt")
out_dir <- "~/scRNA/B/10HC11LN/DEG/data/"
obj <- DGEList(counts=cts, group=clinical_data$Group, samples = clinical_data$sample)
obj$samples
head(obj$counts)
#filter out lowly expressed genes
keep <- filterByExpr(obj, group = clinical_data$Group,
             min.count = 30, min.total.count = 300, large.n = 4, min.prop = 0.6)
table(keep)
obj <- obj[keep, , keep.lib.sizes=FALSE]

## Normalisation for RNA composition using TMM (trimmed mean of M-values)
obj <- calcNormFactors(obj)
obj$samples

## Set up the design matrix
Sample <- factor(clinical_data$sample)
IE <- factor(clinical_data$Group)
IE
## HC HC HC HC HC HC HC HC HC HC LN LN LN LN LN LN LN LN LN LN LN

design <- model.matrix(~IE)

## Estimate dispersion (estimates common dispersion, trended dispersions and tagwise dispersions in one run)
obj <- estimateDisp(obj, design)
plotBCV(obj)

## Calculate Differential Expression
# exact test (only for single-factor experiments)
et <- exactTest(obj)
topTags(et)

# Number of up/downregulated genes at 5% FDR
summary(decideTests(et))
       LN-HC
Down    3549
NotSig  6129
Up      4182
write.csv(as.data.frame(topTags(et, n=Inf)), file=paste0(out_dir, "LNvsHC_EdgeR_DEG_immune.combined.csv"))
## Volcano Plot
in.path <- "/public/home/zhaohuanhuan/scRNA/B/10HC11LN/DEG/data/"
#Read EdgeR data
#edger = read.csv(paste0(in.path, "pseudobulk/Tcell_TIG3vsTIG2_EdgeR_samplesums_exactT.csv"))
edger = read.csv(paste0(in.path, "LNvsHC_EdgeR_samplesums_immune.combined.csv"))


edger$X <- as.character(edger$X)
#remove all genes with logCPM < 1.5
edger <- filter(edger, logCPM > 1.5)

####自定义highlight基因######
highlight <- c("FEZ1","AMKAR","CAMK2N1","MYOM2",
"CD40LG","ID3","CD28","CD86","KLRB1","CD80","PDCD1",
"IL1R2","VSIG4","IFI27","CCR2","C2","C3AR1",
"IL13RA1","TLR2","TLR4","S100A8","PTX3","MPO",#
"C1QA","CD163","MAFB",
"BIRC5","FLT3","OLFM4","DEFA3","DEFA4","LTF","ELL2","ADAMTS2"
)

edger$color <- as.character(ifelse(edger$FDR < 0.1 & edger$logFC > 0.5, "#F8766D", 
                      ifelse(edger$FDR < 0.1 & edger$logFC < -0.5, "#00BFC4", "grey")))

#Plot with ggplot
pdf("immune.combined_edgeR_pseudoBulk_sample_DEG_ViolinPlot_Final-font5-log10Pvalue.pdf",height=5,width=5)
p1 <- ggplot(edger) +
  geom_point(aes(logFC, -log10(PValue)), color=edger$color)+
  geom_label_repel(data = subset(edger, X %in% highlight), aes(log2FC, -log10(PValue), label = X), min.segment.length = 0.1,size=5,label.size = 0.08,max.overlaps = Inf)+
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey20")+
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey20")+
  geom_hline(yintercept = -log(0.1), linetype = "dotted", color = "grey20")+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill="NA"))
p1
dev.off()
#####################全部细胞的基因表达量Boxplots#########################################
```{r}
# load average counts
in_dir = "~/scRNA/B/10HC11LN/DEG/data/immune.combined_sample_sum_counts.csv"
counts <- read.csv(in_dir, row.names="X")

###  setup data table ###
rawdat = as.data.table(counts)
#Normalize counts by dividing through library size
rawdat.cpm <- apply(counts,2, function(x) (x/sum(x))*1000000)
tdat = t(rawdat.cpm)
trnames <- row.names(tdat)
tdat <-as.data.table(tdat)
colnames(tdat) = rownames(counts)
tdat[, condition := trnames]
group.list <- c("HC","HC","HC","HC","HC", "HC","HC","HC","HC","HC",
"LN","LN","LN","LN","LN","LN","LN","LN","LN","LN","LN")  
tdat[, group := group.list]
# format the tables
dat = melt(tdat, id.vars=c('condition', 'group'), variable.name='gene', value.name = 'cpm' , variable.factor = FALSE)

# gene list
GOI <- c("BACH2","BATF", "BCL11B", "BCL11A", "PRDM1",  "PRKCB", "SIRPB1", "KLF2", "ZEB2", "KLF10", "IRF1")

# plot together with edger values
edger$gene <- edger$X
edger$FDR_x = paste0('FDR = ',signif(edger$FDR, digits=3))
edger$p = paste0('p = ',signif(edger$PValue, digits=2))
pdf("allcell_DEG_boxplot.pdf",height=10,width=6.5)
p = subset(dat, gene %in% GOI)%>%
  merge(edger, by='gene')%>%
  ggplot(aes(x=gene,y=cpm,color=group))+
  facet_wrap(~gene+FDR_x+p, scales = "free", ncol = 3)+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.75),
  #size = 1
  )+ scale_color_manual(values = c("HC" = "#00BFC4", "LN" = "#F8766D"))+
  #expand_limits(x=0,y=0)+
  theme_bw()+
  theme(axis.line.x = element_line(colour = "black", size = 0.25),
        axis.line.y = element_line(colour = "black", size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill="NA"),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        #strip.background = element_blank(),
        axis.ticks.x=element_blank())
p
dev.off()
###############
# 输入文件
deg <- read.csv(file.path(input_dir,
  "LNvsHC_EdgeR_DEG_immune.combined.csv"),
  row.names = 1,
  check.names = FALSE
)

head(deg)

summary(deg)
colnames(deg)
deg$gene <- rownames(deg)
####不同FDR和logFC阈值统计 DEG数量
deg <- filter(deg, logCPM > 1.5)#########其实不应该用这个去除低表达量基因，画图时才用，这里为了和之前结果保持一致
FDR_cutoff <- c(
  0.05,
  0.01,
  0.001
)

# absolute logFC thresholds
logFC_cutoff <- c(
  0.5,
  1,
  1.5,
  2
)

#############################
# 3. Calculate DEG numbers
#############################

deg_number <- data.frame()


for(fdr in FDR_cutoff){
  
  for(fc in logFC_cutoff){
    
    
    tmp <- edger %>%
      filter(
        FDR < fdr &
        abs(logFC) > fc
      )
    
    
    up_number <- tmp %>%
      filter(logFC > 0) %>%
      nrow()
    
    
    down_number <- tmp %>%
      filter(logFC < 0) %>%
      nrow()
    
    
    deg_number <- rbind(
      deg_number,
      data.frame(
        FDR = fdr,
        logFC_cutoff = fc,
        Direction = "Up",
        Number = up_number
      ),
      
      data.frame(
        FDR = fdr,
        logFC_cutoff = fc,
        Direction = "Down",
        Number = down_number
      )
    )
    
  }
}

deg_number

#############################
# 4. Save table
#############################
write.csv(
  deg_number,
  paste0(
    table_dir,
    "EdgeR_DEG_threshold_sensitivity_number.csv"
  ),
  row.names = FALSE
)