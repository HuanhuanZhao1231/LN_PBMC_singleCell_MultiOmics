###################每种细胞挑选700个##################
immune.combined <- load("~/scRNA/B/10HC11LN/LN_immune.combined.RData")
immune.combined[["predicted.id"]] <- Idents(immune.combined)
levels(factor(immune.combined@meta.data$predicted.id))
#levels(factor(meta$predicted.id))
meta <- immune.combined@meta.data
###提取700个细胞
meta$order<-"aa"
meta[meta$predicted.id %in% "Bcell",]$order <- c(1:dim(meta[meta$predicted.id %in% "Bcell",])[1]);
sample(c(1:dim(meta[meta$predicted.id %in% "Bcell",])[1]),size=700,replace=FALSE)->include;
start1=dim(meta[meta$predicted.id %in% "Bcell",])[1]+1;
end1=start1-1+dim(meta[meta$predicted.id %in% "CD8T",])[1];
meta[meta$predicted.id %in% "CD8T",]$order<-c(start1:end1);
c(include,sample(c(start1:end1),size=700,replace=FALSE))->include;
start2=end1+1;
end2=start2-1+dim(meta[meta$predicted.id %in% "CD4T",])[1];
meta[meta$predicted.id %in% "CD4T",]$order<-c(start2:end2);
c(include,sample(c(start2:end2),size=700,replace=FALSE))->include;
start3=end2+1;
end3=start3-1+dim(meta[meta$predicted.id %in% "NK",])[1];
meta[meta$predicted.id %in% "NK",]$order<-c(start3:end3);
c(include,sample(c(start3:end3),size=700,replace=FALSE))->include;
start4=end3+1;
end4=start4-1+dim(meta[meta$predicted.id %in% "Mono",])[1];
meta[meta$predicted.id %in% "Mono",]$order<-c(start4:end4);
c(include,sample(c(start4:end4),size=700,replace=FALSE))->include;
start5=end4+1;
end5=start5-1+dim(meta[meta$predicted.id %in% "Neu",])[1];
meta[meta$predicted.id %in% "Neu",]$order<-c(start5:end5);
c(include,sample(c(start5:end5),size=700,replace=FALSE))->include;
start6=end5+1;
end6=start6-1+dim(meta[meta$predicted.id %in% "Mega",])[1];
meta[meta$predicted.id %in% "Mega",]$order<-c(start6:end6);
c(include,sample(c(start6:end6),size=700,replace=FALSE))->include;

meta[meta$order %in% include,] -> c
rownames(c) -> selected
subset <- subset(immune.combined,cells = selected)
######################（1）导出特定基因的数据---导出用于CIBERSORT的细胞基因表达矩阵########################
# 读取基因列表文件，导出特定基因的数据
gene_list <- read.table("~/RNAseq/CIBERSORT/markers_19type_ABC.txt",header=FALSE)
selected_genes <- gene_list$V1
####提取子集
ForCibersort <- subset(subset, features = selected_genes)
# 提取表达矩阵
DefaultAssay(ForCibersort) <- "RNA"
#expression_matrix <- GetAssayData(ForCibersort, slot="data")####标准化的数据
expression_matrix_Count <- GetAssayData(ForCibersort, slot="count")

write.table(expression_matrix_Count,sep="\t",quote=F,row.name=T,col.name=Idents(ForCibersort),file="~/count_700cell_marker.txt")
matrix <- read.table("~/RNAseq/CIBERSORT/count_11type_700cell_marker.txt", header = TRUE, sep = "\t", row.names = 1)
# 检查每一列是否全为0
all_zero_cols <- which(apply(matrix, 2, function(x) all(x == 0)))

# 输出结果
if(length(all_zero_cols) > 0) {
  cat("存在全为0的列：", all_zero_cols, "\n")
} else {
  cat("不存在全为0的列\n")
}
matrix <- matrix[, -all_zero_cols]
write.table(matrix,file="~/RNAseq/CIBERSORT/count_11type_700cell_marker.txt", sep ="\t", row.names =TRUE,col.names =TRUE, quote =FALSE)


counts <- read.table("~/RNAseq/alignOut/deseq_out/countsallToBulkCounts.txt",header=TRUE,row.names=1)

write.table(counts,file="~/RNAseq/CIBERSORT/countsall_symbol.txt", sep ="\t", row.names =TRUE,col.names =TRUE, quote =FALSE)


#########################检查是否有全为0的列###########################
matrix <- read.table("~/RNAseq/CIBERSORT/count5_700cell.txt",header=TRUE,row.names=1)

# 检查每一列是否全为0
all_zero_cols <- which(apply(matrix, 2, function(x) all(x == 0)))

# 输出结果
if(length(all_zero_cols) > 0) {
  cat("存在全为0的列：", all_zero_cols, "\n")
} else {
  cat("不存在全为0的列\n")
}
matrix <- matrix[, -all_zero_cols]
write.table(matrix,file="~/RNAseq/CIBERSORT/count5_700cell.txt", sep ="\t", row.names =TRUE,col.names =TRUE, quote =FALSE)


