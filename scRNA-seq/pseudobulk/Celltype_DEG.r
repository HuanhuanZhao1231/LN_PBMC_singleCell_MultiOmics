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

#####对每种celltypefor循环计算DEG
# List of unique cell types
cell_types <- unique(immune.combined@meta.data$celltype_ABC)
output_file <- paste0(out_dir, "celltype_DEG_output_log.txt")
# 启动文件输出
sink(output_file, append = TRUE)  # append = TRUE 用于追加内容到文件
# Create a list to store subsets
immune_subsets <- list()
# Loop through each cell type and create a subset
for(cell_type in cell_types) {
  immune_subsets[[cell_type]] <- subset(immune.combined, celltype_ABC == cell_type)
}
for(cell_type in cell_types) {
       print(cell_type)
  # Extract raw counts for the current cell type
  raw.counts <- immune_subsets[[cell_type]]@assays$RNA@counts
  # Create clinical data (you may need to adapt this part depending on your dataset)
  meta_data <- immune_subsets[[cell_type]]@meta.data
  barcodes <- rownames(meta_data)
  samples <- meta_data$orig.ident
# 将 orig.ident 和 barcode 连接起来生成新的列名
new_colnames <- paste(samples, barcodes, sep = "_")
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
raw.sums <- raw.counts[, (ncol(raw.counts) - 20):ncol(raw.counts)]
write.csv(raw.sums, file = paste0(out.path, cell_type, "_sample_sum_counts.csv") ,row.names = TRUE)
cts <- raw.sums
clinical_data <- read.table("~/scRNA/B/10HC11LN/DEG/data/clinical_data.txt",header=TRUE)
obj <- DGEList(counts=cts, group=clinical_data$Group, samples = clinical_data$sample)
obj$samples
head(obj$counts)
  # Filter out lowly expressed genes
  keep <- filterByExpr(obj, group = clinical_data$Group,
                       min.count = 30, min.total.count = 300, large.n = 4, min.prop = 0.6)
  obj <- obj[keep, , keep.lib.sizes=FALSE]
  # Normalization (TMM)
  obj <- calcNormFactors(obj)
  # Design matrix
  ## Set up the design matrix
Sample <- factor(clinical_data$sample)
IE <- factor(clinical_data$Group)
IE
design <- model.matrix(~IE)
  # Estimate dispersion
  obj <- estimateDisp(obj, design)
  # Differential expression analysis
  et <- exactTest(obj)
print(topTags(et))
# Number of up/downregulated genes at 5% FDR
print(summary(decideTests(et)))
  # Export results
  write.csv(as.data.frame(topTags(et, n=Inf)), file = paste0(out_dir, "LNvsHC_EdgeR_DEG_", cell_type ,".csv"))
}
sink()###终止定向输出到log文件

