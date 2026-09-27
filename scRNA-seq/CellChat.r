library(CellChat)
library(Seurat)
library(ggplot2)
load("~/scRNA/B/10HC11LN/3score_labled15.8immune.combinedCelltypeGroup.RData")
setwd("~/scRNA/B/10HC11LN/cell-cell-commu")
allcolour_RNA <- readRDS("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts/allcolour_RNA.rds")
immune.combined$cellumap <- ifelse(as.character(immune.combined$celltype_ABC) == "ABC", 
                                   as.character(immune.combined$celltype_ABC), 
                                   as.character(immune.combined$celltype))
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
options(future.globals.maxSize = 8000 * 1024^2)
future::plan("multisession", workers = 3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#将推断的蜂窝通信网络提取为数据帧
df.net <- subsetCommunication(cellchat)
##在信号通路水平推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("cellchat_HC.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
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
# 定义自定义的细胞排列顺序
custom_order <- c("CD4 NC","CD4 ET","CD8 NC","CD8 ET","Prolif","NK","NKR","B IN","B Mem","ABC","Plasma","Mono C","Mono NC-I","Mono NC","Neu","LDG","cDC","pDC")
# 修改 CellChat 对象的细胞顺序
#cellchatHL@idents <- factor(cellchatHL@idents, levels = custom_order)
cellchat_H <- netAnalysis_computeCentrality(cellchat_H, slot.name = "netP")
cellchat_L <- netAnalysis_computeCentrality(cellchat_L, slot.name = "netP")
HL.list <- list(HC=cellchat_H,LN=cellchat_L)
cellchatHL <- mergeCellChat(HL.list,add.names=names(HL.list),cell.prefix=TRUE)
#33333333333333333可视化33333333333333
gg1 <- compareInteractions(cellchatHL,show.legend=F,group=c(1,2),measure="count")
gg2 <- compareInteractions(cellchatHL,show.legend=F,group=c(1,2),measure="weight")
p <- gg1+gg2
ggsave("overview_number_strength.pdf",p,width=6,height=4)
#数量与强度差异网络图
pdf("cellchat_HLgroup_compare.pdf",width=15,height=5)
weight.max <- getMaxWeight(HL.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(HL.list)) {
  netVisual_circle(HL.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(HL.list)[i]))
}
dev.off()
###对比主要的发出者和接受者
#先对每个单独的分组进行计算netAnalysis_computeCentrality
#cellchat_H <- netAnalysis_computeCentrality(cellchat_H, slot.name = "netP")
#cellchat_M <- netAnalysis_computeCentrality(cellchat_M, slot.name = "netP")
#cellchat_S <- netAnalysis_computeCentrality(cellchat_S, slot.name = "netP")
#再运行下面
num.link <- sapply(HL.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(HL.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(HL.list[[i]], title = names(HL.list)[i], weight.MinMax = weight.MinMax,color.use = allcolour_RNA)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("cellchat_HLgroup_source_target.pdf",width=10,height=5)
patchwork::wrap_plots(plots = gg)
dev.off()
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("cellchat_HLgroup_source_target——diff.pdf",width=4.5,height=5)
netAnalysis_diff_signalingRole_scatter(HL.list, comparison = c(1, 2),color.use=allcolour_RNA)
dev.off()
###################两组分组比较#############################
#HM.list <- list(Moderate=cellchat_M,HC=cellchat_H)
#HS.list <- list(Severe=cellchat_S,HC=cellchat_H)
#MS.list <- list(Severe=cellchat_S,Moderate=cellchat_M)
HL.list <- list(HC=cellchat_H,LN=cellchat_L)
#cellchatHM <- mergeCellChat(HM.list,add.names=names(HM.list),cell.prefix=TRUE)
#cellchatHS <- mergeCellChat(HS.list,add.names=names(HS.list),cell.prefix=TRUE)
#cellchatMS <- mergeCellChat(MS.list,add.names=names(MS.list),cell.prefix=TRUE)
cellchatHL <- mergeCellChat(HL.list,add.names=names(HL.list),cell.prefix=TRUE)
##相互作用数量和强度的统计
pdf("cellchat_HLgroup_sum.pdf",width=10,height=5)
gg1 <- compareInteractions(cellchatHL, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchatHL, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
###Differential number of interactions or interaction strength among different cell populations
pdf("cellchat_HLgroup_LNtoHC.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatHL, weight.scale = T)
netVisual_diffInteraction(cellchatHL, weight.scale = T, measure = "weight")
dev.off()
#热图
pdf("cellchat_HLgroup_heatmap.pdf",width=10,height=5)
gg1 <- netVisual_heatmap(cellchatHL)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatHL, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()
#
#pdf("cellchat_HLgroup_scatter_macro.pdf",width=15,height=4)
#gg1 <- netAnalysis_signalingChanges_scatter(cellchatHL, idents.use = "Macro", signaling.exclude = "MIF")
#gg2 <- netAnalysis_signalingChanges_scatter(cellchatHL, idents.use = "Macro", signaling.exclude = c("MIF"))
#patchwork::wrap_plots(plots = list(gg1,gg2))
#dev.off()
###对比每条信号途径的总体信息流###list中数据集在前面的是红色
pdf("cellchat_HLgroup_flow2.pdf",width=10,height=10)
gg1 <- rankNet(cellchatHL, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatHL, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
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
pdf("cellchat_HLcompare_ABC_outgoing.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = 10, targets.use = c(1:19),  comparison = c(1, 2), angle.x = 45)
dev.off()
###
pdf("cellchat_HLcompare_CD8ET_incoming.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 4,  comparison = c(1, 2), angle.x = 45)
dev.off()
pdf("cellchat_HLcompare_CD8ET_incoming_2facet.pdf",height=9,width=15)
gg1 <- netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LN", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LN", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()
pdf("cellchat_HLcompare_Prolif_incoming.pdf",height=9,width=7)
netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 5,  comparison = c(1, 2), angle.x = 45)
dev.off()
pdf("cellchat_HLcompare_Prolif_incoming_2facet.pdf",height=9,width=15)
gg1 <- netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LN", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchatHL, sources.use = c(1:19), targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LN", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

