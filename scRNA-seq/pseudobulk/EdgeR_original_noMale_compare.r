
############################################################
## Female_only sensitivity analysis
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
plot_dir <- "~/scRNA/B/10HC11LN/DEG/Female_only/Results/Figures"
table_dir <- "~/scRNA/B/10HC11LN/DEG/Female_only/Results/Tables"
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
  "hc3",
  "hc5"
)

cts_female <- cts[, !colnames(cts) %in% remove_samples]
clinical_data_female <- clinical_data[!clinical_data$sample %in% remove_samples,]

obj <- DGEList(counts=cts_female, group=clinical_data_female$Group, samples = clinical_data_female$sample)
obj$samples
head(obj$counts)
#filter out lowly expressed genes
keep <- filterByExpr(obj, group = clinical_data_female$Group,
             min.count = 30, min.total.count = 300, large.n = 4, min.prop = 0.6)
table(keep)
FALSE  TRUE
 7577 13906
obj <- obj[keep, , keep.lib.sizes=FALSE]

## Normalisation for RNA composition using TMM (trimmed mean of M-values)
obj <- calcNormFactors(obj)
obj$samples

## Set up the design matrix
Sample <- factor(clinical_data_female$sample)
IE <- factor(clinical_data_female$Group)
IE
###
[1] HC HC HC HC HC HC HC HC LN LN LN LN LN LN LN LN LN LN LN
Levels: HC LN

design <- model.matrix(~IE)

## Estimate dispersion (estimates common dispersion, trended dispersions and tagwise dispersions in one run)
obj <- estimateDisp(obj, design)
pdf(file.path(plot_dir,"edgeR_no_maleHC_plotBCV.pdf"))
plotBCV(obj)
dev.off()
## Calculate Differential Expression
# exact test (only for single-factor experiments)
et <- exactTest(obj)
topTags(et)

# Number of up/downregulated genes at 5% FDR
> summary(decideTests(et))
       LN-HC
Down    3174
NotSig  6958
Up      3774
write.csv(as.data.frame(topTags(et, n=Inf)), file=paste0(out_dir, "LNvsHC_EdgeR_DEG_immune.combined_no_male.csv"))

#################################################
################
deg0 <- read.csv(file.path(input_dir,
  "LNvsHC_EdgeR_DEG_immune.combined.csv"),
  row.names=1,
  check.names = FALSE
)
deg0$gene <-row.names(deg0)
#> dim(deg0)
[1] 13860     5

#################################################
## DEG summary
#################################################
deg_female <- read.csv(file.path(input_dir,
  "LNvsHC_EdgeR_DEG_immune.combined_no_male.csv"),
  row.names = 1,
  check.names = FALSE
)
deg_female$gene <- rownames(deg_female)

> dim(deg_female)
[1] 13906     5

############Pvalue correlation plot########
plot_df <-deg0 %>%dplyr::select(gene,PValue)%>%dplyr::rename(final=PValue)%>%
dplyr::left_join(deg_female %>%dplyr::select(gene,PValue)%>%
dplyr::rename(female=PValue),by="gene")
cor.test(
-log10(plot_df$final),
-log10(plot_df$female)
)

data:  -log10(plot_df$final) and -log10(plot_df$female)
t = 844.56, df = 13834, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9901191 0.9907533
sample estimates:
      cor
0.9904414
###########
ggplot(
plot_df,
aes(
x=-log10(final),
y=-log10(female)
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
ggsave(file.path(plot_dir,"Original_noMale-correlation_log10pvalue.pdf"))
############logFC correlation plot########
plot_df2 <-
deg0 %>%

dplyr::select(
gene,
logFC
)%>%

dplyr::rename(
final=logFC
)%>%

dplyr::left_join(

deg_female %>%

dplyr::select(
gene,
logFC
)%>%

dplyr::rename(
female=logFC
),

by="gene"

)


cor.test(
plot_df2$final,
plot_df2$female
)

data:  plot_df2$final and plot_df2$female
t = 1325.4, df = 13834, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9959529 0.9962134
sample estimates:
      cor
0.9960853
###########################
ggplot(
plot_df2,
aes(
x=final,
y=female
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
y="Female_only log2FC"
)
ggsave(file.path(plot_dir,"Original_Nomale-correlation_FC.pdf"))

####不同FDR和logFC阈值统计 DEG数量
summary_table <- data.frame()


for(fdr in c(0.05,0.01)){

for(lfc in c(0.5,1,1.5,2)){


up <- sum(
deg_female$FDR < fdr &
deg_female$logFC > lfc,
na.rm=TRUE
)


down <- sum(
deg_female$FDR < fdr &
deg_female$logFC < -lfc,
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
"No_male_DEG_summary.csv"
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
deg_female_DEG <- get_DEG(deg_female)

write.csv(
deg_female_DEG,
file.path(table_dir,"res_female_only_edgeR_DEGlist_sameMethod.csv"),
row.names=FALSE
)


deg_final <- get_DEG(deg0)
> length(intersect(deg_female_DEG$gene,deg_final$gene))
[1] 3015

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

gsea_female <-
run_GSEA(deg_female_DEG)

write.csv(
gsea_original2,
file.path(table_dir,"edgeR_DEG_GSEA_original_result.csv"),
row.names=FALSE
)
write.csv(
gsea_female,
file.path(table_dir,"edgeR_methodGSEA_female_only_result.csv"),
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
gsea_female,
"female"
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

ggsave(file.path(plot_dir,"First_method_original_Female_only_GSEA_compare.pdf"),width=7,height=4)

#######

write.csv(
gsea_compare,
file.path(table_dir,"First_method_GSEA_female_only_sensitivity.csv"),
row.names=FALSE
)
library(DESeq2)
#1. 根据gender生成新的三分类变量
cts <- read.csv("~/scRNA/B/10HC11LN/DEG/data/immune.combined_sample_sum_counts.csv", row.names="X")
metadata <- read.table("~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/clinical_data_Covariates.txt",header=T)
rownames(metadata) <- metadata$sample
metadata$gender_group <- NA
metadata$gender_group[
    metadata$Group == "HC"&
    metadata$Gender == 0
] <- "HC_female"
metadata$gender_group[
    metadata$Group == "HC" &
    metadata$Gender == 1
] <- "HC_male"
metadata$gender_group[
    metadata$Group == "LN" &
    metadata$Gender == 0
] <- "LN_female"

metadata$gender_group <- factor(
    metadata$gender_group,
    levels=c(
        "HC_female",
        "HC_male",
        "LN_female"
    )
)
table(metadata$gender_group)


##########
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = metadata,
    design = ~ gender_group
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
    intgroup="gender_group",
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
        color=gender_group,
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


ggsave("~/scRNA/B/10HC11LN/DEG/Female_only/Results/Figures/Sex_group_PCA.pdf",p)
#####################################Venn图可视化#####
library(dplyr)
library(VennDiagram)
library(grid)

# 添加方向
deg_final <- deg_final %>%
  mutate(direction = ifelse(logFC > 0, "UP", "Down"))

deg_female_DEG <- deg_female_DEG %>%
  mutate(direction = ifelse(logFC > 0, "UP", "Down"))


# 提取gene列表
deg_final_UP <- deg_final %>%
  filter(direction == "UP") %>%
  pull(gene)

deg_final_Down <- deg_final %>%
  filter(direction == "Down") %>%
  pull(gene)


deg_female_UP <- deg_female_DEG %>%
  filter(direction == "UP") %>%
  pull(gene)

deg_female_Down <- deg_female_DEG %>%
  filter(direction == "Down") %>%
  pull(gene)

deg_number <- data.frame(
  Dataset = c("deg_final", "deg_final",
              "deg_female_DEG", "deg_female_DEG"),
  Direction = c("UP", "Down",
                "UP", "Down"),
  Number = c(
    length(deg_final_UP),
    length(deg_final_Down),
    length(deg_female_UP),
    length(deg_female_Down)
  )
)

deg_number
> deg_number
         Dataset Direction Number
1      deg_final        UP   1939
2      deg_final      Down   1285
3 deg_female_DEG        UP   1866
4 deg_female_DEG      Down   1243

pdf(file.path(plot_dir,"Original_noMale_Venn_UP.pdf"), width = 5, height = 5)

venn_UP <- venn.diagram(
  x = list(
    deg_final_UP = deg_final_UP,
    deg_female_UP = deg_female_UP
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


pdf(file.path(plot_dir,"Original_noMale_Venn_DOWN.pdf"), width = 5, height = 5)

venn_Down <- venn.diagram(
  x = list(
    deg_final_Down = deg_final_Down,
    deg_female_Down = deg_female_Down
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
    length(setdiff(deg_final_UP, deg_female_UP)),
    length(setdiff(deg_final_Down, deg_female_Down))
  ),
  deg_female_specific = c(
    length(setdiff(deg_female_UP, deg_final_UP)),
    length(setdiff(deg_female_Down, deg_final_Down))
  ),
  Shared = c(
    length(intersect(deg_final_UP, deg_female_UP)),
    length(intersect(deg_final_Down, deg_female_Down))
  )
)

overlap_summary

overlap_summary
  Direction deg_final_specific deg_female_specific Shared
1        UP                126                  53   1813
2      Down                 83                  41   1202

#####################
##########################################################
deg_final_UP_specific <- setdiff(
  deg_final_UP,
  deg_female_UP
)

length(deg_final_UP_specific)

head(deg_final_UP_specific)
deg_final_Down_specific <- setdiff(
  deg_final_Down,
  deg_female_Down
)

length(deg_final_Down_specific)

head(deg_final_Down_specific)

write.table(
  deg_final_UP_specific,
  file = "deg_final_specific_UP_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


write.table(
  deg_final_Down_specific,
  file = "deg_final_specific_Down_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)