
############################################################
## Medication sensitivity analysis
## DESeq2 + apeglm
##因为前面尝试将Steroid——high和dose加入计算Deseq2，但是共线性过高，只能去除Steroid——high的样本
############################################################
library(edgeR)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(pheatmap)
############################################################
## 1. Path setting
############################################################
input_dir <- 
"~/scRNA/B/10HC11LN/DEG/data"
plot_dir <- "~/scRNA/B/10HC11LN/DEG/Medication/Results/Figures"
table_dir <- "~/scRNA/B/10HC11LN/DEG/Medication/Results/Tables"
############################################################
## 2. Clinical information
############################################################
####
cts <- read.csv("~/scRNA/B/10HC11LN/DEG/data/immune.combined_sample_sum_counts.csv", row.names="X")
clinical_data <- read.table("~/scRNA/B/10HC11LN/DEG/data/clinical_data.txt",header=TRUE)
rownames(clinical_data) <- clinical_data$sample
out_dir <- "~/scRNA/B/10HC11LN/DEG/data/"
cts <- cts[
,
clinical_data$sample
]


stopifnot(
all(
colnames(cts)==rownames(clinical_data)
)
)


remove_samples <- c(
  "pbmc1k",
  "pbmc4k",
  "pbmc8k",
  "pbmc11k"
)

cts_low <- cts[, !colnames(cts) %in% remove_samples]
clinical_data_low <- clinical_data[!clinical_data$sample %in% remove_samples,]

obj <- DGEList(counts=cts_low, group=clinical_data_low$Group, samples = clinical_data_low$sample)
obj$samples
head(obj$counts)
#filter out lowly expressed genes
keep <- filterByExpr(obj, group = clinical_data_low$Group,
             min.count = 30, min.total.count = 300, large.n = 4, min.prop = 0.6)
table(keep)
FALSE  TRUE
 7564 13919
obj <- obj[keep, , keep.lib.sizes=FALSE]

## Normalisation for RNA composition using TMM (trimmed mean of M-values)
obj <- calcNormFactors(obj)
obj$samples

## Set up the design matrix
Sample <- factor(clinical_data_low$sample)
IE <- factor(clinical_data_low$Group)
IE
## HC HC HC HC HC HC HC HC HC HC LN LN LN LN LN LN LN LN LN LN LN

design <- model.matrix(~IE)

## Estimate dispersion (estimates common dispersion, trended dispersions and tagwise dispersions in one run)
obj <- estimateDisp(obj, design)

plotBCV(obj)
ggsave(file.path(plot_dir,"edgeR_no_high_steroid_plotBCV.pdf"))
## Calculate Differential Expression
# exact test (only for single-factor experiments)
et <- exactTest(obj)
topTags(et)

# Number of up/downregulated genes at 5% FDR
> summary(decideTests(et))
       LN-HC
Down    3219
NotSig  6893
Up      3807
write.csv(as.data.frame(topTags(et, n=Inf)), file=paste0(out_dir, "LNvsHC_EdgeR_DEG_immune.combined_no_Steroid_high.csv"))

#################################################


#################################################
## DEG summary
#################################################
deg_low <- read.csv(file.path(input_dir,
  "LNvsHC_EdgeR_DEG_immune.combined_no_Steroid_high.csv"),
  row.names = 1,
  check.names = FALSE
)
deg_low$gene <- rownames(deg_low)
############logFC correlation plot
plot_df <-deg0 %>%dplyr::select(gene,PValue)%>%rename(final=PValue)%>%
left_join(deg_low %>%dplyr::select(gene,PValue)%>%
rename(low=PValue),by="gene")
cor.test(
-log10(plot_df$final),
-log10(plot_df$low)
)

data:  -log10(plot_df$final) and -log10(plot_df$low)
t = 435.5, df = 13820, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9642936 0.9665587
sample estimates:
      cor
0.9654443
###########
ggplot(
plot_df,
aes(
x=-log10(final),
y=-log10(low)
)
)+

geom_point(
alpha=0.3,
size=0.5
)+

geom_smooth(
method="lm"
)+

theme_classic()+

labs(
x="Original -log10(Pvalue)",
y="Steroid-low log2FC"
)
ggsave(file.path(plot_dir,"Original_steroid-low-correlation_log10pvalue.pdf"))
################################
plot_df2 <-
deg0 %>%

dplyr::select(
gene,
logFC
)%>%

rename(
final=logFC
)%>%

left_join(

deg_low %>%

dplyr::select(
gene,
logFC
)%>%

rename(
low=logFC
),

by="gene"

)


cor.test(
plot_df2$final,
plot_df2$low
)

data:  plot_df2$final and plot_df2$low
t = 893.06, df = 13820, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9911583 0.9917264
sample estimates:
     cor
0.991447
###########################
ggplot(
plot_df2,
aes(
x=final,
y=low
)
)+

geom_point(
alpha=0.3,
size=0.5
)+

geom_smooth(
method="lm"
)+

theme_classic()+

labs(
x="Original log2FC",
y="Steroid-low log2FC"
)
ggsave(file.path(plot_dir,"Original_steroid-low-correlation_FC.pdf"))

####不同FDR和logFC阈值统计 DEG数量

summary_table <- data.frame()


for(fdr in c(0.05,0.01)){

for(lfc in c(0.5,1,1.5,2)){


up <- sum(
deg_low$FDR < fdr &
deg_low$logFC > lfc,
na.rm=TRUE
)


down <- sum(
deg_low$FDR < fdr &
deg_low$logFC < -lfc,
na.rm=TRUE
)


summary_table <-
rbind(
summary_table,
data.frame(
FDR=fdr,
LFC=lfc,
UP=up,
DOWN=down
)
)


}

}


write.csv(

summary_table,

file.path(
table_dir,
"No_steroid_high_DEG_summary.csv"
),

row.names=FALSE

)
############
get_DEG <- function(x){
x %>%
filter(
FDR<0.05,
abs(logFC)>0.5
)

}

deg_low_DEG <- get_DEG(deg_low)

write.csv(
deg_low_DEG,
file.path(table_dir,"res_low_steroid_edgeR_DEGlist_sameMethod.csv"),
row.names=FALSE
)
################
deg0 <- read.csv(file.path(input_dir,
  "LNvsHC_EdgeR_DEG_immune.combined.csv"),
  row.names=1,
  check.names = FALSE
)
deg0$gene <-row.names(deg0)

deg_final <- get_DEG(deg0)
write.csv(
deg_final,
file.path(table_dir,"res_final_DEGlist.csv"),
row.names=FALSE
)
> length(intersect(deg_low_DEG$gene,deg_final$gene))
[1] 2956


###############GSEA富集分析########
library(msigdbr)
library(clusterProfiler)


hallmark <- msigdbr(
species="Homo sapiens",
category="H"
)


TERM2GENE <- hallmark[,c(
"gs_name",
"gene_symbol"
)]

run_GSEA <- function(res){
geneList <- res$logFC

names(geneList)<-
rownames(res)

geneList <-
sort(
geneList,
decreasing=TRUE
)

GSEA(
geneList,
TERM2GENE=TERM2GENE
)

}

gsea_original2 <-
run_GSEA(deg_final)

gsea_low_steroid2 <-
run_GSEA(deg_low_DEG)

write.csv(
gsea_original2,
file.path(table_dir,"edgeR_DEG_GSEA_original_result.csv"),
row.names=FALSE
)
write.csv(
gsea_low_steroid2,
file.path(table_dir,"edgeR_methodGSEA_low_steroid_result.csv"),
row.names=FALSE
)
###############12. 比较IFN/TLR/T cell pathway
pathway_interest <- c(

"HALLMARK_INTERFERON_ALPHA_RESPONSE",

"HALLMARK_INTERFERON_GAMMA_RESPONSE",

"HALLMARK_TNFA_SIGNALING_VIA_NFKB",

"HALLMARK_T_CELL_RECEPTOR_SIGNALING",

"HALLMARK_IL6_JAK_STAT3_SIGNALING"

)



extract_GSEA <- function(x,name){


as.data.frame(x) %>%

filter(
ID %in% pathway_interest
)%>%

dplyr::select(
ID,
NES,
p.adjust
)%>%

mutate(
model=name
)

}



gsea_compare <- bind_rows(

extract_GSEA(
gsea_original2,
"Original"
),

extract_GSEA(
gsea_low_steroid2,
"Steroid low"
)

)

gsea_compare

ggplot(
gsea_compare,
aes(
x=model,
y=NES,
fill=ID
)
)+
geom_col(
position="dodge"
)+
theme_classic()

ggsave(file.path(plot_dir,"First_method_GSEA_compare.pdf"),width=7,height=4)

#######

write.csv(
gsea_compare,
file.path(table_dir,"First_method_GSEA_medication_sensitivity.csv"),
row.names=FALSE
)
###################看Steroid_high对转录的影响——PCA画图#####
library(DESeq2)
#1. 根据Steroid_high生成新的三分类变量
rownames(metadata) <- metadata$sample
metadata$Steroid_group <- NA
metadata$Steroid_group[
    metadata$disease == "HC"
] <- "HC"
metadata$Steroid_group[
    metadata$disease == "LN" &
    metadata$Steroid_high == 0
] <- "LN_low"
metadata$Steroid_group[
    metadata$disease == "LN" &
    metadata$Steroid_high == 1
] <- "LN_high"
metadata$Steroid_group <- factor(
    metadata$Steroid_group,
    levels=c(
        "HC",
        "LN_low",
        "LN_high"
    )
)
table(metadata$Steroid_group)


##########
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ Steroid_group
)


dds <- estimateSizeFactors(dds)


vsd <- vst(
    dds,
    blind = TRUE
)


vst_mat <- assay(vsd)
##########PCA计算#######
library(ggplot2)


pcaData <- plotPCA(
    vsd,
    intgroup="Steroid_group",
    returnData=TRUE
)


percentVar <- round(
    100 * attr(pcaData,"percentVar")
)


head(pcaData)
###########绘图#
p <- ggplot(
    pcaData,
    aes(
        PC1,
        PC2,
        color=Steroid_group,
        label=name
    )
)+
geom_point(
    size=4
)+
geom_text_repel(
    size=3
)+
xlab(
    paste0(
        "PC1: ",
        percentVar[1],
        "%"
    )
)+
ylab(
    paste0(
        "PC2: ",
        percentVar[2],
        "%"
    )
)+
theme_classic()+
theme(
    legend.title=element_blank()
)


ggsave("~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/Medication/Results/Figures/Steroid_group_PCA.pdf",p)
#####################################Venn图可视化#####
library(dplyr)
library(VennDiagram)
library(grid)

# 添加方向
deg_final <- deg_final %>%
  mutate(direction = ifelse(logFC > 0, "UP", "Down"))

deg_low_DEG <- deg_low_DEG %>%
  mutate(direction = ifelse(logFC > 0, "UP", "Down"))


# 提取gene列表
deg_final_UP <- deg_final %>%
  filter(direction == "UP") %>%
  pull(gene)

deg_final_Down <- deg_final %>%
  filter(direction == "Down") %>%
  pull(gene)


deg_low_UP <- deg_low_DEG %>%
  filter(direction == "UP") %>%
  pull(gene)

deg_low_Down <- deg_low_DEG %>%
  filter(direction == "Down") %>%
  pull(gene)

deg_number <- data.frame(
  Dataset = c("deg_final", "deg_final",
              "deg_low_DEG", "deg_low_DEG"),
  Direction = c("UP", "Down",
                "UP", "Down"),
  Number = c(
    length(deg_final_UP),
    length(deg_final_Down),
    length(deg_low_UP),
    length(deg_low_Down)
  )
)

deg_number
> deg_number
      Dataset Direction Number
1   deg_final        UP   1939
2   deg_final      Down   1285
3 deg_low_DEG        UP   1881
4 deg_low_DEG      Down   1297

pdf("Venn_UP.pdf", width = 5, height = 5)

venn_UP <- venn.diagram(
  x = list(
    deg_final_UP = deg_final_UP,
    deg_low_UP = deg_low_UP
  ),
  filename = NULL,
  fill = c("#E64B35", "#4DBBD5"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = 0
)

grid.draw(venn_UP)

dev.off()


pdf("Venn_DOWN.pdf", width = 5, height = 5)

venn_Down <- venn.diagram(
  x = list(
    deg_final_Down = deg_final_Down,
    deg_low_Down = deg_low_Down
  ),
  filename = NULL,
  fill = c("#00A087", "#3C5488"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = 0
)


grid.draw(venn_Down)

dev.off()

overlap_summary <- data.frame(
  Direction = c("UP", "Down"),
  deg_final_specific = c(
    length(setdiff(deg_final_UP, deg_low_UP)),
    length(setdiff(deg_final_Down, deg_low_Down))
  ),
  deg_low_specific = c(
    length(setdiff(deg_low_UP, deg_final_UP)),
    length(setdiff(deg_low_Down, deg_final_Down))
  ),
  Shared = c(
    length(intersect(deg_final_UP, deg_low_UP)),
    length(intersect(deg_final_Down, deg_low_Down))
  )
)

overlap_summary

> overlap_summary
  Direction deg_final_specific deg_low_specific Shared
1        UP                150               92   1789
2      Down                118              130   1167