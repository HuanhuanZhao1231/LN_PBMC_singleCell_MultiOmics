library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)
library(DESeq2)
#library(apeglm)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("~/scRNA/B/10HC11LN")
load("immune.combined.RData")
meta <- immune.combined@meta.data
#导出细分细胞亚群的计数#####
composition_count <- meta %>%
group_by(
orig.ident,
celltype_ABC
)%>%
summarise(
n=n(),
.groups="drop"
)
#转宽：
composition_matrix <- composition_count %>%
pivot_wider(
names_from = celltype_ABC,
values_from=n,
values_fill=0
)

composition_matrix <- as.data.frame(composition_matrix)

rownames(composition_matrix)<-
composition_matrix$orig.ident


composition_matrix$orig.ident <- NULL


composition_matrix

out_dir <- "~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/"
write.csv(
  composition_matrix,
  file = paste0(out_dir, "immune_cell_composition_matrix_Celltype_ABC.csv"),
  row.names = TRUE
)


#Step 1. 生成major cell type（用于scCODA）

#不要直接用19个cell type。

#原因：

#pDC=227个细胞

#会影响模型。

#定义：
meta$major_celltype <- NA


meta$major_celltype[
meta$celltype_ABC %in%
c("CD4ET","CD4NC","CD8ET","CD8NC","Prolif")
] <- "T"


meta$major_celltype[
meta$celltype_ABC %in%
c("ABC","BIN","BMem","Plasma")
] <- "B"


meta$major_celltype[
meta$celltype_ABC %in%
c("NK","NKR")
] <- "NK"


meta$major_celltype[
meta$celltype_ABC %in%
c(
"MonoC",
"MonoNC",
"MonoNCI",
"Neu",
"LDG",
"cDC",
"pDC")
] <- "Myeloid"


meta$major_celltype[
meta$celltype_ABC=="Mega"
] <- "Mega"


table(meta$major_celltype)#########

###########################
###生成major cell type（用于scCODA）
meta$major_celltype <- NA


meta$major_celltype[
meta$celltype_ABC %in%
c("CD4ET","CD4NC")
] <- "CD4T"

meta$major_celltype[
meta$celltype_ABC %in%
c("CD8ET","CD8NC","Prolif")
] <- "CD8T"

meta$major_celltype[
meta$celltype_ABC %in%
c("ABC","BIN","BMem","Plasma")
] <- "B"


meta$major_celltype[
meta$celltype_ABC %in%
c("NK","NKR")
] <- "NK"


meta$major_celltype[
meta$celltype_ABC %in%
c(
"MonoC",
"MonoNC",
"MonoNCI",
"cDC",
"pDC"
)###DC数量太少，不适合单列
] <- "Mono"


meta$major_celltype[
meta$celltype_ABC %in%
c(
"Neu",
"LDG"
)
] <- "Neu"

meta$major_celltype[
meta$celltype_ABC=="Mega"
] <- "Mega"


table(meta$major_celltype)


#Step 2. 生成cell composition matrix

#这是scCODA输入。

#每个sample每种细胞数量
composition_count <- meta %>%
group_by(
orig.ident,
major_celltype
)%>%
summarise(
n=n(),
.groups="drop"
)
#转宽：
composition_matrix <- composition_count %>%
pivot_wider(
names_from = major_celltype,
values_from=n,
values_fill=0
)

composition_matrix <- as.data.frame(composition_matrix)

rownames(composition_matrix)<-
composition_matrix$orig.ident


composition_matrix$orig.ident <- NULL


composition_matrix

out_dir <- "~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/"
write.csv(
  composition_matrix,
  file = paste0(out_dir, "immune_cell_composition_matrix_majorCelltype.csv"),
  row.names = TRUE
)

composition_fraction <- 
  composition_matrix / rowSums(composition_matrix)

write.csv(
  composition_fraction,
  file = paste0(out_dir, "immune_cell_composition_fraction_majorCelltype.csv"),
  row.names = TRUE
)
#4. 我建议你顺便保存一个带Group信息的版本

#因为后面 scCODA 和 ggplot 都需要：
clinical <- read.table(
  "~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/clinical_data_Covariates.txt",
  header = TRUE,
  sep = "\t"
)

rownames(clinical) <- clinical$sample


composition_metadata <- cbind(
  clinical[rownames(composition_matrix),],
  composition_matrix
)


write.csv(
  composition_metadata,
  file=paste0(out_dir,
              "immune_cell_composition_with_clinical.csv"),
  row.names=TRUE
)