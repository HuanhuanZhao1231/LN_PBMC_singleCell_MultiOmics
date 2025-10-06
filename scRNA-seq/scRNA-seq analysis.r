
library(Seurat)
library(harmony)
library(dplyr)
setwd("~/scRNA/B/10HC11LN")
pbmc.list<-list(hc2,hc3,hc4,hc5,hc6,hc7,hc8,hc9,hc10,hc11,pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10,pbmc11)
pbmc.list<-lapply(X =pbmc.list, FUN =function(x){
x<-ScaleData(x, features =features, verbose =FALSE)
x<-RunPCA(x, features =features, verbose =FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
save(immune.combined,file="integrated_immune.combined.Rdata")
immune.combined <- ScaleData(immune.combined, verbose = FALSE)###必须先缩放再PCA
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(immune.combined), npcs = 30)
#run harmony
immune.combined <- RunHarmony(immune.combined, group.by.vars = 'orig.ident',project.dim = F)######不加project.dim = F 会报错Error in data.use %*% cell.embeddings : non-conformable arguments
save(immune.combined,file="afterHarmony_immune.combined.Rdata")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
immune.combined <- JackStraw(immune.combined, num.replicate = 100)
immune.combined <- ScoreJackStraw(immune.combined, dims = 1:20)
pdf("elbowplot.pdf")
ElbowPlot(immune.combined)
dev.off()
load("afterHarmony_immune.combined.Rdata")
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
for (i in c(6,7,8,9,10,11,12,13,14)) {
    for (k in c(20,40,60)) {
immune.combined <-immune.combined %>% 
RunUMAP(reduction = "harmony", n.neighbors=22,min.dist = 0.2, dims = 1:i) %>% 
FindNeighbors(reduction = "harmony", dims = 1:i,k.param =k) ####默认k.param =20，根据需要可改为其他
immune.combined <- FindClusters(immune.combined,resolution =0.8,graph.name = "RNA_snn")###不加graph.name = "RNA_snn"报错 Provided graph.name not present in Seurat object
# Visualization
pdfname <- paste0("integrationmin0.2","dim",i,"k",k,".pdf")
p1 <- DimPlot(immune.combined, reduction = "umap",group.by="orig.ident",raster=FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
p <- p1 + p2
ggsave(pdfname,plot=p,width=13,height=5)
}
    }
immune.combined <-immune.combined %>% 
RunUMAP(reduction = "harmony", n.neighbors=22,min.dist = 0.2, dims = 1:15) %>% 
FindNeighbors(reduction = "harmony", dims = 1:15,k.param =20) ####默认k.param =20，根据需要可改为其他
immune.combined <- FindClusters(immune.combined,resolution = 0.8,graph.name = "RNA_snn")###不加graph.name = "RNA_snn"报错 Provided graph.name not present in Seurat object
# Visualization
pdf("integration_15.8K20umap.pdf",height=5,width=13)
p1 <- DimPlot(immune.combined, reduction = "umap",group.by="orig.ident",raster=FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
p <- p1+p2
p
dev.off()
#########IFNScore-Inflammation response score##
##-------------------AUCell-------------------------
####准备矩阵###
library(GSEABase)
library(AUCell)
library(DelayedArray)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
exprMatrix <- immune.combined@assays$RNA@data
###########准备GeneSets##########
gene_list <- read.table("~/scRNA/B/integration/5sample/data/ISG_Nature_immunity_heatmap_100gene.txt",header=FALSE)
isg_genes <- gene_list$V1
#gene_list2 <- read.table("~/scRNA/B/integration/5sample/data/Inflammatory.txt",header=FALSE)
#infla_genes <- gene_list3$V1
geneSets <- GeneSet(unique(isg_genes), setName="geneSet1")
geneSets
cells_rankings <- AUCell_buildRankings(exprMatrix)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.05)#
AUCell_auc <- as.numeric(getAUC(cells_AUC))
immune.combined@meta.data[['ISGscore']] <- AUCell_auc
####添加分组信息
group <- immune.combined@meta.data$orig.ident
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c("pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k','pbmc5k','pbmc3k','pbmc9k','pbmc2k')){
    group[which(group==i_group)] <- 'LN'}
unique(group)##可以看group的值
immune.combined@meta.data[['Group']] <- group
#AIgroup
group <- immune.combined@meta.data$orig.ident
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c("pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k')){
    group[which(group==i_group)] <- 'Moderate'
}
for(i_group in c('pbmc5k','pbmc3k','pbmc9k','pbmc2k')){
    group[which(group==i_group)] <- 'Severe'
}####按顺序生成一个group列表
unique(group)##可以看group的值
immune.combined@meta.data[['AI2Group']] <- group
#############################分样本画#######################
pdf("ISGscore_AUCell_AI2Group_p.pdf",width=4,height=5)
data <- FetchData(immune.combined,vars = c("ISGscore","AI2Group"))
data$orig.ident <- factor(data$orig.ident, levels = custom_order, ordered = TRUE)
data <- data[order(data$orig.ident), ]
my_comparisons <- list( c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Moderate", "Severe"))
#可以查看计算的p值stat.test <- data %>% t_test(ISGscore ~ AI2Group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
p <- ggplot(data = data,aes(AI2Group,ISGscore,fill = AI2Group))
p + geom_boxplot(outlier.size=0.2) + theme_bw() + RotatedAxis() + labs(title = "ISGscore",y = "Score") + 
theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 10,face = "bold"),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12)) +
coord_cartesian(ylim = c(0, 0.4)) +
stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.signif",label.y = c(0.34,0.36,0.38))
#############cell-cell communication#################
library(CellChat)
library(Seurat)
library(ggplot2)
load("~/scRNA/B/10HC11LN/3score_labled15.8immune.combinedCelltypeGroup.RData")
setwd("~/scRNA/B/10HC11LN/cell-cell-commu")
allcolour_RNA <- readRDS("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts/allcolour_RNA.rds")
#2.2 分析HC样本的细胞通讯网络
HC.combined <- subset(immune.combined,subset = Group == "Healthy")
data.HC.input = HC.combined@assays$RNA@data # normalized data matrix
HC.meta = HC.combined@meta.data
unique(HC.meta$cellumap)
cellchat <- createCellChat(object = data.HC.input, meta = HC.meta, group.by = "cellumap")
###Add cell information into meta slot of the object (Optional)
cellchat <- addMeta(cellchat, meta = HC.meta)
cellchat <- setIdent(cellchat, ident.use = "cellumap") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
####Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
options(future.globals.maxSize = 8000 * 1024^2)##不加这个会报错
future::plan("multisession", workers = 3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#将推断的蜂窝通信网络提取为数据帧
df.net <- subsetCommunication(cellchat)#返回一个数据帧，其中包含配体/受体水平上所有推断的细胞间通讯。设置slot.name = "netP"为访问信号通路级别的推断通信
##在信号通路水平推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#我们还可以可视化聚合的细胞间通信网络。例如，使用圆形图显示任意两个细胞组之间的相互作用次数或总相互作用强度（权重）。
groupSize <- as.numeric(table(cellchat@idents))
pdf("cellchat_HC.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#由于细胞间的通信网络复杂，我们可以检查每个细胞组发送的信令。这里我们还控制参数，edge.weight.max以便我们可以比较不同网络之间的边权重。
mat <- cellchat@net$weight
pdf("cellchat_HC_cell.pdf",height=20,width=20)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_HC.rds")
#2.3 分析LN样本的细胞通讯网络
LN.combined <- subset(immune.combined,subset = Group == "LN")
data.LN.input = LN.combined@assays$RNA@data # normalized data matrix
LN.meta = LN.combined@meta.data
unique(LN.meta$cellumap)
cellchat <- createCellChat(object = data.LN.input, meta = LN.meta, group.by = "cellumap")
###Add cell information into meta slot of the object (Optional)
cellchat <- addMeta(cellchat, meta = LN.meta)
cellchat <- setIdent(cellchat, ident.use = "cellumap") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
####Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#将推断的蜂窝通信网络提取为数据帧
df.net <- subsetCommunication(cellchat)#返回一个数据帧，其中包含配体/受体水平上所有推断的细胞间通讯。设置slot.name = "netP"为访问信号通路级别的推断通信
##在信号通路水平推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#我们还可以可视化聚合的细胞间通信网络。例如，使用圆形图显示任意两个细胞组之间的相互作用次数或总相互作用强度（权重）。
groupSize <- as.numeric(table(cellchat@idents))
pdf("cellchat_LN.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#由于细胞间的通信网络复杂，我们可以检查每个细胞组发送的信令。这里我们还控制参数，edge.weight.max以便我们可以比较不同网络之间的边权重。
mat <- cellchat@net$weight
pdf("cellchat_LN_cell.pdf",height=20,width=20)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, file = "cellchat_LN.rds")
#2.4 合并cellchat对象
cellchat_H <- readRDS("cellchat_HC.rds")
cellchat_L <- readRDS("cellchat_LN.rds")
##对比每种细胞的传入传出信号
library(ComplexHeatmap)
i = 1
pdf("cellchat_HL_outgoing_pattern.pdf",width=15,height=16)
pathway.union <- union(HL.list[[i]]@netP$pathways, HL.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(HL.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(HL.list)[i], width = 10, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(HL.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(HL.list)[i+1], width = 10, height = 15)
ht1+ht2
dev.off()
pdf("cellchat_HL_incoming_pattern.pdf",width=15,height=16)
pathway.union <- union(HL.list[[i]]@netP$pathways, HL.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(HL.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(HL.list)[i], width = 10, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(HL.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(HL.list)[i+1], width = 10, height = 15)
ht1+ht2
dev.off()
# 定义自定义的细胞排列顺序
custom_order <- c("CD4 NC","CD4 ET","CD8 NC","CD8 ET","Prolif","NK","NKR","B IN","B Mem","ABC","Plasma","Mono C","Mono INT","Mono NC","Neu1","Neu2","cDC","pDC")
###
pdf("cellchat_HLcompare_Prolif_outgoing.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = 5, targets.use = c(1:19),  comparison = c(1, 2), angle.x = 45)
dev.off()
pdf("cellchat_HLcompare_cDC_outgoing.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = 17, targets.use = c(1:19),  comparison = c(1, 2), angle.x = 45)
dev.off()
pdf("cellchat_HLcompare_Plasma_outgoing.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = 11, targets.use = c(1:19),  comparison = c(1, 2), angle.x = 45)
dev.off()
###
pdf("cellchat_HLcompare_CD8ET_incoming.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 4,  comparison = c(1, 2), angle.x = 45)
dev.off()
#########