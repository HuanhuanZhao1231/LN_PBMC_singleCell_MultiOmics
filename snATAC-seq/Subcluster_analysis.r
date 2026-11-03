####以myeloid cells为例###Bcell和T/NK也进行亚群细分#######
#8.1将Seurat对象转换为SingleCellExperiment对象
load("~/scRNA/B/10HC11LN/LN_immune.combinedCelltypeGroup.RData")
Lym <- subset(immune.combined,idents=c("CD4 NC","CD4 ET","CD8 NC","CD8 ET","NK","NKR","Prolif"))
Mye <- subset(immune.combined,idents=c("Mono C","Mono NC","Mono NC-I","Neu","LDG","pDC","cDC"))
DefaultAssay(Lym) <- "RNA"
celltype <- Idents(Lym)
Lym@meta.data[['celltype']] <- celltype
Lym.sce <- as.SingleCellExperiment(Lym)
colnames(colData(Lym.sce))
table(colData(Lym.sce)$celltype)
saveRDS(Lym.sce, file="~/snATAC/B/ArchR/LN_7type_Lym.sce.rds")
#Mye亚群分析##
Mye <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_ProjHeme5_Myeloid")
Mye <- addIterativeLSI(
    ArchRProj = Mye,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    clusterParams = list(resolution=c(2), sampleCells=30000, 
    #maxClusters=6, 
    n.start=10),
    sampleCellsPre = 30000,
    varFeatures = 25000,
    dimsToUse = 1:25,
    force = TRUE
)
Mye <- addHarmony(
    ArchRProj = Mye,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
Mye <- addClusters(
    input = Mye,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
)
set.seed(1)
Mye <- addUMAP(
    ArchRProj = Mye, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = Mye, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Mye, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Mye_Plot-UMAP-Sample_Cluster_0.4_35_0.4_test.pdf", ArchRProj = Mye, addDOC = FALSE, width = 10, height = 10)
saveArchRProject(Mye)
# Make various cluster plots:
Mye <- addImputeWeights(Mye)
Mye_markerGenes <- c(
"CD14","LYZ","FCGR3A",    
"FCER1A","CST3","CD68",
"LILRA4","IRF7","TCF4","PLD4",#pDC
"CD1C","ITGAX",#cDC
"RETN","FCGR2A","S100A8","IL1B","CSF3R","CMTM2",#Neu
"FCGR3B","FUT4"#LDG
)
p2 <- plotEmbedding(
    ArchRProj = Mye, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = Mye_markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(Mye)
)
p2c <- lapply(p2, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5,legendTextSize = 3) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
plotPDF(plotList = p2c, 
    name = "Mye_Plot-UMAP-Marker-Genes-RNA0.4_0.4-W-Imputation_integration_genescore.pdf", 
    ArchRProj = Mye, 
    addDOC = FALSE, width = 5, height = 5)
###去除C1(n<100,离群)
idxPass <- which(Mye$Clusters %in% c("C2","C3","C4","C5","C6","C7"))
cellsPass <- Mye$cellNames[idxPass]
Mye_rm1 <- Mye[cellsPass, ]
Mye_rm1 <- addIterativeLSI(
    ArchRProj = Mye_rm1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    clusterParams = list(resolution=c(2), sampleCells=30000, 
    #maxClusters=6, 
    n.start=10),
    sampleCellsPre = 30000,
    varFeatures = 25000,
    dimsToUse = 1:25,
    force = TRUE
)
Mye_rm1 <- addHarmony(
    ArchRProj = Mye_rm1,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
Mye_rm1 <- addClusters(
    input = Mye_rm1,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
)
set.seed(1)
Mye_rm1 <- addUMAP(
    ArchRProj = Mye_rm1, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = Mye_rm1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Mye_rm1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Myerm1_Plot-UMAP-Sample_Cluster_0.4_35_0.4.pdf", ArchRProj = Mye_rm1, addDOC = FALSE, width = 10, height = 10)
saveArchRProject(Mye_rm1,outputDirectory="~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Mye_rm1_0.4_0.4",load=FALSE)
####2024.04.29####
#Mye_rm1 <- loadArchRProject("/public/home/zhaohuanhuan/snATAC/B/ArchR/Mye_rm12_0.4_0.4")
Mye.sce <- readRDS("~/snATAC/B/ArchR/7type_Mye.sce.rds")
addArchRThreads(threads = 1)
Mye_rm1 <- addGeneIntegrationMatrix(
    ArchRProj = Mye_rm1, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = Mye.sce, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    force = TRUE,
    groupRNA = "celltype", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "FineClust_RNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore" #Name of column where prediction score from scRNA
)

cM <- as.matrix(confusionMatrix(Mye_rm1$Clusters, Mye_rm1$FineClust_RNA))
####重新命名
cM <- confusionMatrix(Mye_rm1$Clusters, Mye_rm1$predictedGroup)
labelOld <- rownames(cM)
labelOld
#[1] "C5" "C6" "C4" "C2" "C3" "C1"
remapClust <- c(
    "C1" = "Neu",
    "C2" = "DC",
    "C3" = "Mono NC",
    "C4" = "Mono C",
    "C5" = "Mono C",
    "C6" = "Mono NC-I"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelOld, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
#  [1] "GMP"        "B"          "PreB"       "CD4.N"      "Mono"      
#  [6] "Erythroid"  "Progenitor" "CD4.M"      "pDC"        "NK"        
# [11] "CLP"        "Mono"
#[1] "Bmem" "Bmem" "BIN"  "PC"   "ABC"  "Bmem"
Mye_rm1$Clusters2 <- mapLabels(Mye_rm1$Clusters, newLabels = labelNew2, oldLabels = labelOld)
pdf(file=paste0(archrdir,"/Plots/labeledumap.pdf"),height=5,width=5)
p2 <- plotEmbedding(ArchRProj = Mye_rm1, colorBy = "cellColData", pal = allcolour_atac ,name = "Clusters2", embedding = "UMAP")
p2
dev.off()
###用定义的monocyte细胞亚群的名称替换大群的未细分的mono的名称
# 获取proj的行数
proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_ProjHeme4")
num_rows <- nrow(getCellColData(proj))###
# 遍历proj的每一行
for (i in 1:num_rows) {
  # 获取proj的当前行的rowname
  current_rowname <- rownames(getCellColData(proj))[i]  
  # 检查当前行的rowname是否在Mye_rm1中存在
  if (current_rowname %in% rownames(getCellColData(Mye_rm1))) {
    # 获取Mye_rm1中与当前行rowname相匹配的行的索引
    matching_row_index <- match(current_rowname, rownames(getCellColData(Mye_rm1)))    
    # 获取对应行的内容
    matching_row_content <- getCellColData(Mye_rm1)[matching_row_index, ]    
    # 将proj$FineClust的内容替换为matching_row_content
    proj$FineClust[i] <- matching_row_content$Clusters2
  }
}
###
for (i in 1:nrow(getCellColData(proj))) {
if (proj$FineClust[i] == "Mono") {
    proj$FineClust[i] <- "Neu"
  }
}
pdf(file=paste0(archrdir,"/Plots/projHeme5_MonoNC_labeledumap.pdf"),height=5,width=5)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", pal = allcolour_atac ,name = "FineClust", embedding = "UMAP")
p2
dev.off()
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-proj_N9_Monosub_ProjHeme4", load = TRUE)




























#peak可视化
p <- plotBrowserTrack(
    ArchRProj = Mye_rm12, 
    groupBy = "Clusters2", 
    geneSymbol = B_markerGenes, 
    upstream = 50000,
    downstream = 50000
)
plotPDF(plotList = p, 
    name = "Mye_rm12_Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = Mye_rm12, 
    addDOC = FALSE, width = 5, height = 5)
#### Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
###添加无约束积分#耗时久100min左右，建议提交运行
#addArchRThreads(threads = 1)##默认不是1个线程，会内存不足
B.sce <- readRDS("/public/home/zhaohuanhuan/snATAC/B/ArchR/B.sce_from_PBMC.B.combined_labled_umap13r2.rds")
Mye_rm12 <- addGeneIntegrationMatrix(
    ArchRProj = Mye_rm12, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",###2024.04.30好像要和前面一致用harmony
    seRNA = B.sce,
    addToArrow = FALSE,
    sampleCellsATAC = 30000,
    sampleCellsRNA = 30000,
    groupRNA = "celltype2",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)
#From scRNA
cBIN <- "TrB|BIN|BIN2|BIN3|ISGhigh-BIN"
cBMem <- "Bmem|Bmem2|Bmem3|ABC|ISGhigh-Bmem"
cP <- "PB|PC"
#Assign scATAC to these categories
clustBIN <- c("C5")
clustBIN
clustBmem <- c("C2","C3","C4","C6")
clustBmem
clustP <- c("C1")
clustP
##
rnaBIN <- colnames(B.sce)[grep(cBIN, colData(B.sce)$celltype2)]
rnaBmem <- colnames(B.sce)[grep(cBMem, colData(B.sce)$celltype2)]
rnaP <- colnames(B.sce)[grep(cP, colData(B.sce)$celltype2)]
head(rnaBIN)
groupList <- SimpleList(
    BIN = SimpleList(
        ATAC = Mye_rm12$cellNames[Mye_rm12$Clusters %in% clustBIN],
        RNA = rnaBIN
    ),
    BMem = SimpleList(
        ATAC = Mye_rm12$cellNames[Mye_rm12$Clusters %in% clustBmem],
        RNA = rnaBmem
    ),
    P = SimpleList(
        ATAC = Mye_rm12$cellNames[Mye_rm12$Clusters %in% clustP],
        RNA = rnaP
    )   
)

#We pass this list to the `groupList` parameter of the `addGeneIntegrationMatrix()` function to constrain our integration. Note that, in this case, we are still not adding these results to the Arrow files (`addToArrow = FALSE`). We recommend checking the results of the integration thoroughly against your expectations prior to saving the results in the Arrow files. We illustrate this process in the next section of the book.
#~30 minutes
####每个 scATAC-seq 细胞添加Pseudo-scRNA-seq profiles
Mye_rm12 <- addGeneIntegrationMatrix(
    ArchRProj = Mye_rm12, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",###2024.04.30好像要和前面一致用harmony
    seRNA = B.sce,
    addToArrow = TRUE,
    sampleCellsATAC = 30000,
    sampleCellsRNA = 30000,
    force= TRUE,
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
saveArchRProject(Mye_rm12, load = TRUE)
#Making Pseudo-bulk Replicates
Mye_rm12 <- addGroupCoverages(ArchRProj = Mye_rm12, groupBy = "Clusters2")
#Our scRNA labels
table(Mye_rm12$Clusters2)
pathToMacs2 <- "/public/home/zhaohuanhuan/.conda/envs/macs2/bin/macs2"
Mye_rm12 <- addReproduciblePeakSet(
    ArchRProj = Mye_rm12, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(Mye_rm12)
projHemeTmp <- addReproduciblePeakSet(
    ArchRProj = Mye_rm12, 
    groupBy = "Clusters2",
    peakMethod = "Tiles",
    method = "p"
)
getPeakSet(projHemeTmp)
saveArchRProject(ArchRProj = Mye_rm12, outputDirectory = "Save-Mye_rm12_addPeak", load = FALSE)
projHeme5 <- addPeakMatrix(Mye_rm12)
getAvailableMatrices(projHeme5)
###Identifying Marker Peaks with ArchR
#Our scRNA labels
table(projHeme5$Clusters2)
markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
markerList$ABC
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
markerList$ABC
# 原始顺序后面做聚类，没有排序，根据第一次的结果进行调整，调整后的列名顺序
new_order <- c("PC","BIN", "Bmem", "ABC")
# 根据new_order重新排序列名
new_colnames <- colnames(markersPeaks)[order(match(colnames(markersPeaks), new_order))]
# 根据新的列名顺序重新设置SummarizedExperiment对象的列
markersPeaks <- markersPeaks[, new_colnames]
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)###transpose决定图是横的还是竖的
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Mye_rm12_Peak-Marker-Heatmap_sort2", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
###marker_Peak MA 和火山图
pma <- markerPlot(seMarker = markersPeaks, name = "PC", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
pv <- markerPlot(seMarker = markersPeaks, name = "PC", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "PC-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)