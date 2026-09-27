library(Seurat)
library(ggplot2)
load("~/scRNA/B/10HC11LN/immune.combinedGroup_celltype_ABC_cellumap.RData")
setwd("~/scRNA/B/10HC11LN/QC")
####################画TSS和log10nFrag的小提琴图######################
# 1. 
cellData <- as.data.frame(immune.combined@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA","percent.mt")])
# 3. 小提琴图
p1 <- ggplot(cellData, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 14) +
  labs(y = "nCount_RNA", x = "orig.ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 ggsave(plot=p1,"scRNA-seq_21_QC_nCount_RNA.pdf",width=13,height=5)
# 4. 画 nFeature_RNA 小提琴图
p2 <- ggplot(cellData, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 14) +
  labs(y = "nFeature_RNA", x = "orig.ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot=p2,"scRNA-seq_21_QC_nFeature_RNA.pdf",width=13,height=5)
# 4. 画 percent.mt 小提琴图
p3 <- ggplot(cellData, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 14) +
  labs(y = "percent.mt", x = "orig.ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot=p3,"scRNA-seq_21_QC_percent.mt.pdf",width=13,height=5)