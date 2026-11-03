suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(ggrastr)
})
addArchRGenome("hg19")
subproj <- loadArchRProject("./filtered_output")
############
subProj <- addIterativeLSI(
    ArchRProj = subProj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 3,
    sampleCellsFinal =50000,
    projectCellsPre = TRUE, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = 1, 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 1:15,
       force = TRUE #之前运行过降维，不加force,会报错已存在
)
###用harmony降维###
subProj <- addHarmony(
    ArchRProj = subProj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
subProj <- addClusters(
    input = subProj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
)###这里resolution决定分群数量
subProj <- addUMAP(
    ArchRProj = subProj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.3, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = subProj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_LSI3dim15_r1R0.4n35_minD0.4.pdf", ArchRProj = subProj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = subProj, outputDirectory = "Save-subProj_11_LSI3dim15", load = FALSE)
##########
head(subProj$Clusters)
table(subProj$Clusters)
cM <- confusionMatrix(paste0(subProj$Clusters), paste0(subProj$Sample))
cM
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
p
plotPDF(p, name = "cluster_pheatmap_0.2_35_0.3.pdf", ArchRProj = subProj, addDOC = FALSE, width = 5, height = 5)
# Identifying Marker Genes
markersGS <- getMarkerFeatures(
    ArchRProj = subProj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
markerGenes_all  <- c("CD34","SOX4",#Progenitor
"STMN1","TOP2A",#Prolif
"CD3E","CD4","CD40LG","IL7R","TNFRSF4","RTKN2","FOXP3",
"CD8A",#
"CCR7",
"PRF1","GZMH","GZMK",
"NCAM1","GZMA","GZMB","KLRB1","GNLY","NKG7","FCGR3A",
"MS4A1","CD79A",
"CR2","CXCR5",#
"CD38","CD9",#TrB
"TNFRSF17","PCNA","MKI67",
"IFIT1","RSAD2","ISG15","OAS3",
"TBX21","ITGAX","FCRL5","ZEB2","CD19","CD27",
"SDC1","CD38",
"NCAM1", "NKG7","KLRB1",
"CD19", "CD14", "FCGR3A","LYZ","CST3",
"FCGR3A","FCGR3B","CEBPA","CSF3R","CMTM2", "S100A8", "RETN","ITGAM", 
"FCGR2A","FCER1A","CST3","CD68","LILRA4","CD1C","HBB","PPBP")
markerGenes_heatmap  <- c("CD34","SOX4",#Progenitor
"STMN1","TOP2A",#Prolif
"CD3E","CD4","CD40LG","IL7R",
"CD8A",
"CCR7","GZMH",
"NCAM1","KLRB1",
"MS4A1","CD19","CD27",
"CD38","TNFRSF17","PCNA","MKI67",
"IFIT1","OAS3",
"SDC1",
"CD14", "FCGR3A","LYZ","CST3",
"FCGR3A","FCGR3B", "S100A8", "RETN","ITGAM", 
"FCER1A","CD68","LILRA4","CD1C","HBB","PPBP")
###
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes_heatmap,
  transpose = FALSE
) 
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap_heatmap_test", width = 8, height = 6, ArchRProj = subProj, addDOC = FALSE)
#subProj <- addGeneScoreMatrix(subProj, force=TRUE)
##7.4Visualizing Marker Genes on an Embedding
p <- plotEmbedding(
    ArchRProj = subProj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes_all, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = subProj, 
    addDOC = FALSE, width = 5, height = 5)
##使用MAGIC填充标记基因,根据邻近细胞填充基因得分对信号进行平滑化处理
subProj <- addImputeWeights(subProj)
p <- plotEmbedding(
    ArchRProj = subProj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes_all, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(subProj)
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation_R0.4.pdf", 
    ArchRProj = subProj, 
    addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = subProj, outputDirectory = "Save-subProj_11_LSI3_60_0.4", load = FALSE)
addArchRThreads(threads = 1)##
imm.sce <- readRDS("~/snATAC/B/ArchR/LN_18type_imm.sce.rds")
##无拘束整合#########
subProj <- addGeneIntegrationMatrix(
    ArchRProj = subProj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = imm.sce,
    addToArrow = FALSE,
    sampleCellsATAC = 30000,
    sampleCellsRNA = 30000,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)
saveArchRProject(subProj)
##查看无拘束整合分群信息######
cM <- as.matrix(confusionMatrix(subProj$Clusters, subProj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
#From scRNA
cTNK <- "CD4|CD8|NK"
cMye <- "Mono|Neu|LDG|Mega"
cB <- "B|Plasma"
#Assign scATAC to these categories
clustMye <- c("C1","C2","C3","C4","C5","C6","C7")
clustMye
clustTNK <- c("C8","C9","C10","C11","C12","C13","C14","C15","C16")
clustTNK
clustB <- c("C17","C18","C19","C20")
clustB
##
rnaTNK <- colnames(imm.sce)[grep(cTNK, colData(imm.sce)$celltype)]
rnaMye <- colnames(imm.sce)[grep(cMye, colData(imm.sce)$celltype)]
rnaB <- colnames(imm.sce)[grep(cB, colData(imm.sce)$celltype)]
head(rnaTNK)
groupList <- SimpleList(
    TNK = SimpleList(
        ATAC = subProj$cellNames[subProj$Clusters %in% clustTNK],
        RNA = rnaTNK
    ),
    B = SimpleList(
        ATAC = subProj$cellNames[subProj$Clusters %in% clustB],
        RNA = rnaB
    ),
    Mye = SimpleList(
        ATAC = subProj$cellNames[subProj$Clusters %in% clustMye],
        RNA = rnaMye
    )   
)

#We pass this list to the `groupList` parameter of the `addGeneIntegrationMatrix()` function to constrain our integration. Note that, in this case, we are still not adding these results to the Arrow files (`addToArrow = FALSE`). We recommend checking the results of the integration thoroughly against your expectations prior to saving the results in the Arrow files. We illustrate this process in the next section of the book.
#~30 minutes
#subProj <- loadArchRProject("/public/home/zhaohuanhuan/snATAC/B/ArchR/Save-subProj_11_LSI3dim15")
subProj <- addGeneIntegrationMatrix(
    ArchRProj = subProj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = imm.sce,
    addToArrow = FALSE, 
    sampleCellsATAC = 30000,
    sampleCellsRNA = 30000,
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)
#比较无约束和约束积分
pal <- paletteDiscrete(values = colData(imm.sce)$celltype)
pal
p1 <- plotEmbedding(
    subProj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)
p2 <- plotEmbedding(
    subProj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal
)
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration_test.pdf", ArchRProj = subProj, addDOC = FALSE, width = 5, height = 5)
#保存
saveArchRProject(ArchRProj = subProj, outputDirectory = "Save-subProj_11_LSI3dim15_immsce", load = FALSE)
####2024.04.23####
####每个 scATAC-seq 细胞添加Pseudo-scRNA-seq profiles
#3h
projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = subProj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = imm.sce,
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
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "/public/home/zhaohuanhuan/snATAC/B/ArchR/Save-ProjHeme3_2_11_LSI3dim15", load = FALSE)
######
projHeme3 <- addUMAP(
    ArchRProj = projHeme3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 60, #增加会使图更紧密
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ClustersLSI3_landmark5w_0.3_60_0.4_beforeLable.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 10, height = 10)
####重新命名
cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

remapClust <- c(
    "C1" = "aMy4",
    "C2" = "aMy1",
    "C3" = "aMy2",
    "C4" = "aMy3",
    "C14" = "aMy5",
    "C15" = "aMy6",
    "C10" = "aBc1",
    "C11" = "aBc2",
    "C5" = "aTc1",
    "C6" = "aTc2",
    "C7" = "aTc3",
    "C8" = "aTc4",
    "C9" = "aTc5",
    "C12" = "aTc6",
    "C13" = "aTc7"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
pdf(file="~/snATAC/B/ArchR/Save-ProjHeme3_2_11_LSI3dim15/Plots/labeledumap.pdf",height=5,width=5)
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p2
dev.off()
###
subClusterGroups <- list(
  "T/NK" = c("aTc1","aTc2","aTc3","aTc4","aTc5","aTc6","aTc7"), 
  "Myeloid" = c("aMy1","aMy2","aMy3","aMy4","aMy5","aMy6"),
  "Bcells" = c("aBc1","aBc2")
  )
# subgroups are now not non-overlapping
subClusterCells <- lapply(subClusterGroups, function(x){
  getCellNames(projHeme3)[as.character(projHeme3@cellColData[["Clusters2"]]) %in% x]
  })

subClusterArchR <- function(projHeme3, subCells, outdir){
  # Subset an ArchR project for focused analysis

  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
      ArchRProj = projHeme3,
      cells = subCells,
      outputDirectory = outdir,
      dropCells = TRUE,
      force = TRUE
  )
  saveArchRProject(sub_proj)
}

subgroups <- c("T/NK", "Myeloid","Bcells")

# Generate and cluster each of the subprojects
sub_proj_list <- lapply(subgroups, function(sg){
  message(sprintf("Subsetting %s...", sg))
  outdir <- sprintf("~/snATAC/B/ArchR/subclustered_%s", sg)
  subClusterArchR(projHeme3, subCells=subClusterCells[[sg]], outdir=outdir)
})
names(sub_proj_list) <- subgroups
saveArchRProject(ArchRProj = sub_proj_list$T/NK,outputDirectory="~/snATAC/B/ArchR/subclustered_T/NK",load=FALSE)
saveArchRProject(ArchRProj = sub_proj_list$Myeloid,outputDirectory="~/snATAC/B/ArchR/subclustered_Myeloid",load=FALSE)
saveArchRProject(ArchRProj = sub_proj_list$Bcells,outputDirectory="~/snATAC/B/ArchR/subclustered_Bcells",load=FALSE)
# Save project
saveArchRProject(projHeme3)