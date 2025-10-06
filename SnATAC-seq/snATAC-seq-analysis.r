suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(ggrastr)
})
addArchRGenome("hg19")
inputFiles <- list.files(pattern = "\\.gz$")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = Samples,
  minTSS = 0, # Don't filter at this point
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)
#构建Arrow Project
proj <-ArchRProject(
ArrowFiles = ArrowFiles,
outputDirectory ="./unfiltered_output_TSS0",
copyArrows = TRUE
) 
proj <- loadArchRProject("./unfiltered_output_TSS0")
addArchRThreads(threads = 24)
idxPass <- which(proj$TSSEnrichment >= 8)
cellsPass <- proj$cellNames[idxPass]
proj <- proj[cellsPass,]
#########
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)
png("TSS-vs-Frags.png")
plot(p)
dev.off()
#每个样本
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
###
proj <- addTileMatrix(proj, force=TRUE)
proj <- addGeneScoreMatrix(proj, force=TRUE)
proj <- addDoubletScores(
  proj,
    dimsToUse=1:20, 
    scaleDims=TRUE, 
    LSIMethod=2
)
saveArchRProject(proj)
proj <- filterDoublets(proj, filterRatio = 1)
# Reduce Dimensions with Iterative LSI (<5 minutes)
set.seed(1)
proj <- loadArchRProject("./proj_TSS8_Tile_Gene")
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    sampleCellsPre = 15000,
    varFeatures = 50000, 
    dimsToUse = 1:25,
    force = TRUE
)
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
# Identify Clusters from Iterative LSI
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2,
    force = TRUE
)
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 40, 
    minDist = 0.2, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ClustersLSI2_TSS8_minD0.2_40.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 40, 
    minDist = 0.2, 
    metric = "cosine",
    force = TRUE
)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP-Sample-Clustersharmony_LSI2.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)
# Relabel clusters so they are sorted by cluster size
proj <- relabelClusters(proj)
proj <- addImputeWeights(proj)
# Make various cluster plots:
proj <- visualizeClustering(proj, pointSize=pointSize, sampleCmap=sample_cmap, diseaseCmap=disease_cmap)

# Save filtered ArchR project
saveArchRProject(proj)
##########################################################################################
# Visualize Data
##########################################################################################
# Now, identify likely cells:
identifyCells <- function(df, TSS_cutoff=6, nFrags_cutoff=2000, minTSS=5, minFrags=1000, maxG=4){
    # Identify likely cells based on gaussian mixture modelling.
    # Assumes that cells, chromatin debris, and other contaminants are derived from
    # distinct gaussians in the TSS x log10 nFrags space. Fit a mixture model to each sample
    # and retain only cells that are derived from a population with mean TSS and nFrags passing
    # cutoffs
    ####################################################################
    # df = data.frame of a single sample with columns of log10nFrags and TSSEnrichment
    # TSS_cutoff = the TSS cutoff that the mean of a generating gaussian must exceed
    # nFrags_cutoff = the log10nFrags cutoff that the mean of a generating gaussian must exceed
    # minTSS = a hard cutoff of minimum TSS for keeping cells, regardless of their generating gaussian
    # maxG = maximum number of generating gaussians allowed

    cellLabel <- "cell"
    notCellLabel <- "not_cell"

    if(nFrags_cutoff > 100){
        nFrags_cutoff <- log10(nFrags_cutoff)
        minFrags <- log10(minFrags)
    } 
    
    # Fit model
    set.seed(1)
    mod <- Mclust(df, G=2:maxG, modelNames="VVV")

    # Identify classifications that are likely cells
    means <- mod$parameters$mean

    # Identify the gaussian with the maximum TSS cutoff
    idents <- rep(notCellLabel, ncol(means))
    idents[which.max(means["TSSEnrichment",])] <- cellLabel

    names(idents) <- 1:ncol(means)

    # Now return classifications and uncertainties
    df$classification <- idents[mod$classification]
    df$classification[df$TSSEnrichment < minTSS] <- notCellLabel
    df$classification[df$nFrags < minFrags] <- notCellLabel
    df$cell_uncertainty <- NA
    df$cell_uncertainty[df$classification == cellLabel] <- mod$uncertainty[df$classification == cellLabel]
    return(list(results=df, model=mod))
}

# Run classification on all samples
minTSS <- 5
samples <- unique(proj$Sample)
cellData <- getCellColData(proj)
cellResults <- lapply(samples, function(x){
  df <- cellData[cellData$Sample == x,c("nFrags","TSSEnrichment")]
  df$log10nFrags <- log10(df$nFrags)
  df <- df[,c("log10nFrags","TSSEnrichment")]
  identifyCells(df, minTSS=minTSS)
  })
names(cellResults) <- samples

# Save models for future reference
saveRDS(cellResults, file = paste0(proj@projectMetadata$outputDirectory, "/cellFiltering.rds"))

# Plot filtering results
for(samp in samples){
    df <- as.data.frame(cellResults[[samp]]$results)
    cell_df <- df[df$classification == "cell",]
    non_cell_df <- df[df$classification != "cell",]

    xlims <- c(log10(500), log10(100000))
    ylims <- c(0, 18)
    # QC Fragments by TSS plot w/ filtered cells removed:
    p <- ggPoint(
        x = cell_df[,1], 
        y = cell_df[,2], 
        size = 1.5,
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = xlims,
        ylim = ylims,
        title = sprintf("%s droplets plotted", nrow(cell_df)),
        rastr = TRUE
    )
    # Add grey dots for non-cells
    p <- p + geom_point_rast(data=non_cell_df, aes(x=log10nFrags, y=TSSEnrichment), color="light grey", size=0.5)
    p <- p + geom_hline(yintercept = minTSS, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
    plotPDF(p, name = paste0(samp,"_EM_model_filtered_cells_TSS-vs-Frags.pdf"), ArchRProj = proj, addDOC = FALSE)
}

# Real cells pass QC filter and for C_SD_POOL are classified singlets
realCells <- getCellNames(proj)[(proj$cellCall == "cell")]
#& (proj$DemuxletClassify %ni% c("AMB", "DBL")) & (proj$Sample2 != "C_SD_POOL")]
subProj <- subsetArchRProject(proj, cells=realCells, 
    outputDirectory="./filtered_output", dropCells=TRUE, force=TRUE)

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

# Add Infered Doublet Scores to ArchR project (~5-10 minutes)
#subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)###LN12,报错：LN12 (11 of 12) : Correlation of UMAP Projection is below 0.9 (normally this is ~0.99)
subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2,force=TRUE)
# Visualize numeric metadata per grouping with a violin plot now that we have created an ArchR Project.
plotList <- list()
plotList[[1]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "TSSEnrichment"
)
plotList[[2]] <- plotGroups(ArchRProj = subProj, 
  groupBy = "Sample", 
  colorBy = "colData", 
  name = "DoubletEnrichment"
)
plotPDF(plotList = plotList, name = "TSS-Doublet-Enrichment", width = 4, height = 4,  ArchRProj = subProj, addDOC = FALSE)

# Filter doublets:
subProj <- filterDoublets(subProj, filterRatio = 1)

# Save filtered ArchR project
saveArchRProject(subProj)
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
imm.sce <- readRDS("/public/home/zhaohuanhuan/snATAC/B/ArchR/18type_imm.sce.rds")
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
#[1] "C14" "C3"  "C1"  "C15" "C8"  "C5"  "C10" "C4"  "C11" "C13" "C6"  "C12"
#[13] "C9"  "C2"  "C7"
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
#[1] "Neu"       "Mono C"    "Mono NC-I" "Neu"       "CD4 ET"    "CD8 ET"
#[7] "B IN"      "Mono C"    "Plasma"    "NK"        "NK"        "NK"
#[13] "CD4 NC"    "Mono C"    "CD4 NC"
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
pdf(file="/public/home/zhaohuanhuan/snATAC/B/ArchR/Save-ProjHeme3_2_11_LSI3dim15/Plots/labeledumap.pdf",height=5,width=5)
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p2
dev.off()
###
subClusterGroups <- list(
  "Lymphoid" = c("aTc1","aTc2","aTc3","aTc4","aTc5","aTc6","aTc7"), 
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

subgroups <- c("Lymphoid", "Myeloid","Bcells")

# Generate and cluster each of the subprojects
sub_proj_list <- lapply(subgroups, function(sg){
  message(sprintf("Subsetting %s...", sg))
  outdir <- sprintf("/public/home/zhaohuanhuan/snATAC/B/ArchR/subclustered_%s", sg)
  subClusterArchR(projHeme3, subCells=subClusterCells[[sg]], outdir=outdir)
})###不能把子文件夹保存在/public/home/zhaohuanhuan/snATAC/B/ArchR/subCluster文件夹下面，因为是从subCluster文件夹的ArchRproject进行分割，并新建文件，这个路径会覆盖前面的
names(sub_proj_list) <- subgroups
saveArchRProject(ArchRProj = sub_proj_list$Lymphoid,outputDirectory="/public/home/zhaohuanhuan/snATAC/B/ArchR/subclustered_Lymphoid",load=FALSE)
saveArchRProject(ArchRProj = sub_proj_list$Myeloid,outputDirectory="/public/home/zhaohuanhuan/snATAC/B/ArchR/subclustered_Myeloid",load=FALSE)
#saveArchRProject(ArchRProj = sub_proj_list$Bcells,outputDirectory="/public/home/zhaohuanhuan/snATAC/B/ArchR/subclustered_Bcells",load=FALSE)
# Save project
saveArchRProject(projHeme3)
###B细胞亚群分析##
Bcell <- loadArchRProject("/public/home/zhaohuanhuan/snATAC/B/ArchR/subclustered_Bcells")
Bcell <- addIterativeLSI(
    ArchRProj = Bcell,
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
Bcell <- addHarmony(
    ArchRProj = Bcell,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
Bcell <- addClusters(
    input = Bcell,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
)
set.seed(1)
Bcell <- addUMAP(
    ArchRProj = Bcell, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = Bcell, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Bcell, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Bcell_Plot-UMAP-Sample_Cluster_0.4_35_0.4_test2.pdf", ArchRProj = Bcell, addDOC = FALSE, width = 10, height = 10)
saveArchRProject(Bcell)
# Make various cluster plots:
Bcell <- addImputeWeights(Bcell)
B_markerGenes  <- c(
"CD86",#Activated B cells
"CD27","XBP1","TNFRSF17","MZB1","IRF4",#Plasma cells#TNFRSF17是BCMA
"PRDM1","PCNA","MKI67",#plasmablasts
"CD19","SDC1","CD38","CXCR4","IL6R","BCMA","PTPRC",
"IFIT1","RSAD2","ISG15","OAS3",#ISG-high B
"CD79A","CD79B","MS4A1","TCL1A","FCER2","IL4R",#B-immature,MS4A1是CD20
"CD27","TNFRSF13B",#B-memory
"TBX21","ITGAX","FCRL5","ZEB2",#ABCells
"CR2","CXCR5",
"CD9","CD72","PPP1R14A","FAM129C","NEIL1",#Translational Bcell-ShennanScience中差异表达基因
"NT5E","ENTPD1","FAS","FCRL4"#
)
p2 <- plotEmbedding(
    ArchRProj = Bcell, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = B_markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(Bcell)
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
    name = "Bcellrm12_Plot-UMAP-Marker-Genes-RNA0.4_0.4-W-Imputation_integration_genescore2.pdf", 
    ArchRProj = Bcell, 
    addDOC = FALSE, width = 5, height = 5)
###去除C1，C2小群(n<100,离群)
Bcell <- addIterativeLSI(
    ArchRProj = Bcell,
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
Bcell <- addHarmony(
    ArchRProj = Bcell,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
Bcell <- addClusters(
    input = Bcell,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
)
set.seed(1)
Bcell <- addUMAP(
    ArchRProj = Bcell, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = Bcell, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Bcell, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Bcell_Plot-UMAP-Sample_Cluster_0.4_35_0.4.pdf", ArchRProj = Bcell, addDOC = FALSE, width = 10, height = 10)
saveArchRProject(Bcell,outputDirectory="/public/home/zhaohuanhuan/snATAC/B/ArchR/Bcell_0.4_0.4",load=FALSE)
####2024.04.29####
B.sce <- readRDS("/public/home/zhaohuanhuan/snATAC/B/ArchR/B.sce_from_PBMC.B.combined_labled_umap13r2.rds")
Bcell <- addGeneIntegrationMatrix(
    ArchRProj = Bcell, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = B.sce, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    force = TRUE,
    groupRNA = "celltype2", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
    nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "FineClust_RNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore" #Name of column where prediction score from scRNA
)
cM <- as.matrix(confusionMatrix(Bcell$Clusters, Bcell$FineClust_RNA))
####重新命名
cM <- confusionMatrix(Bcell$Clusters, Bcell$predictedGroup)
labelOld <- rownames(cM)
labelOld
#[1] "C3" "C4" "C5" "C1" "C2" "C6"
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
#[1] "Bmem" "Bmem" "BIN"  "PC"   "ABC"  "Bmem"
Bcell$Clusters2 <- mapLabels(Bcell$Clusters, newLabels = labelNew, oldLabels = labelOld)
pdf(file="./labeledumap.pdf",height=5,width=5)
p2 <- plotEmbedding(ArchRProj = Bcell, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p2
dev.off()
#peak可视化
p <- plotBrowserTrack(
    ArchRProj = Bcell, 
    groupBy = "Clusters2", 
    geneSymbol = B_markerGenes, 
    upstream = 50000,
    downstream = 50000
)
plotPDF(plotList = p, 
    name = "Bcell_Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = Bcell, 
    addDOC = FALSE, width = 5, height = 5)
#### Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
###添加无约束积分#耗时久100min左右，建议提交运行
#addArchRThreads(threads = 1)##默认不是1个线程，会内存不足
B.sce <- readRDS("/public/home/zhaohuanhuan/snATAC/B/ArchR/B.sce_from_PBMC.B.combined_labled_umap13r2.rds")
Bcell <- addGeneIntegrationMatrix(
    ArchRProj = Bcell, 
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
        ATAC = Bcell$cellNames[Bcell$Clusters %in% clustBIN],
        RNA = rnaBIN
    ),
    BMem = SimpleList(
        ATAC = Bcell$cellNames[Bcell$Clusters %in% clustBmem],
        RNA = rnaBmem
    ),
    P = SimpleList(
        ATAC = Bcell$cellNames[Bcell$Clusters %in% clustP],
        RNA = rnaP
    )   
)

#We pass this list to the `groupList` parameter of the `addGeneIntegrationMatrix()` function to constrain our integration. Note that, in this case, we are still not adding these results to the Arrow files (`addToArrow = FALSE`). We recommend checking the results of the integration thoroughly against your expectations prior to saving the results in the Arrow files. We illustrate this process in the next section of the book.
#~30 minutes
####每个 scATAC-seq 细胞添加Pseudo-scRNA-seq profiles
Bcell <- addGeneIntegrationMatrix(
    ArchRProj = Bcell, 
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
saveArchRProject(Bcell, load = TRUE)
#Making Pseudo-bulk Replicates
Bcell <- addGroupCoverages(ArchRProj = Bcell, groupBy = "Clusters2")
#Our scRNA labels
table(Bcell$Clusters2)
pathToMacs2 <- "/public/home/zhaohuanhuan/.conda/envs/macs2/bin/macs2"
Bcell <- addReproduciblePeakSet(
    ArchRProj = Bcell, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(Bcell)
projHemeTmp <- addReproduciblePeakSet(
    ArchRProj = Bcell, 
    groupBy = "Clusters2",
    peakMethod = "Tiles",
    method = "p"
)
getPeakSet(projHemeTmp)
saveArchRProject(ArchRProj = Bcell, outputDirectory = "Save-Bcell_addPeak", load = FALSE)
projHeme5 <- addPeakMatrix(Bcell)
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
plotPDF(heatmapPeaks, name = "Bcell_Peak-Marker-Heatmap_sort2", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
###marker_Peak MA 和火山图
pma <- markerPlot(seMarker = markersPeaks, name = "PC", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
pv <- markerPlot(seMarker = markersPeaks, name = "PC", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "PC-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
###浏览器轨迹中的Marker_Peak
#调整细胞顺序
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "Clusters2", 
    geneSymbol = c("TBX21"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["ABC"],
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$ABC)
plotPDF(p, name = "Plot-Tracks-With-Features_ABC", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
##Pairwise Testing Between Groups
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ABC",
  bgdGroups = "Bmem"
)
pma <- markerPlot(seMarker = markerTest, name = "ABC", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma
pv <- markerPlot(seMarker = markerTest, name = "ABC", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "ABC-vs-Bmem-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
####12.Motif Enrichment in Differential Peaks##上面的usegroup VS bdggroup
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif",force=TRUE)
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black",
        max.overlaps = 20###显示标签的个数
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggUp
motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
motifsDo
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black",
        max.overlaps = 20,###能显示的标签数
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggDo
plotPDF(ggUp, ggDo, name = "ABC-vs-Bmem-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
##Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Bcell_Motifs-Enriched-Marker-Heatmap2", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
##ArchR Enrichment
##1.编码TF结合位点
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
#Annotation ArchR-Hg19-v1.Anno does not exist! Downloading..
#Error in download.file(url = url, destfile = file.path(annoPath, basename(url)),  :
#  cannot open URL 'https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Hg19-v1.Anno'
enrichEncode
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "Bcell_EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
###
#projHeme5 <- loadArchRProject(path ="/public/home/zhaohuanhuan/Save-ProjHeme5_motifEnrich")
###与bulk ATAC
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "ATAC")
enrichATAC <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "ATAC",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
enrichATAC
heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapATAC, name = "Bcell_ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
##22#Codex TFBS
projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "Codex")
enrichCodex <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Codex",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
enrichCodex
heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapCodex, name = "Bcell_Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
## 13 ChromVAR Deviatons Enrichment with ArchR
if("Motif" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}
projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
) ##耗时略长，1h
saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-Bcell_ProjHeme5_motif",load = FALSE)
############################
#projHeme5 <- loadArchRProject(path="Save-ProjHeme5_motif")
options(ggrepel.max.overlaps = Inf)####
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Bcell_Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
motifs <- c("Tcf7","TBX21","EOMES","PAX8","EBF1","IRF1","IRF8","CEBPA","ATF4","PAX5","JUNB","KLF4","NR4A1")
#motifs <- c("SPI1","C-JUN","FOS","KLF2","KLF4","CEBPA","ATF4","NR4A1","JUNB")##单核细胞
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
# markerMotifs
# [1] "z:EOMES_788"           "z:TBX21_780"           "z:TCF7L1_763"
# [4] "z:TCF7L2_762"          "z:TCF7_750"            "z:PAX5_709"
# [7] "z:PAX8_707"            "z:IRF8_633"            "z:IRF1_629"
#[10] "z:CEBPA_155"           "z:ATF4_122"            "z:EBF1_67"
#[13] "z:SREBF1_22"           "deviations:EOMES_788"  "deviations:TBX21_780"
#[16] "deviations:TCF7L1_763" "deviations:TCF7L2_762" "deviations:TCF7_750"
#[19] "deviations:PAX5_709"   "deviations:PAX8_707"   "deviations:IRF8_633"
#[22] "deviations:IRF1_629"   "deviations:CEBPA_155"  "deviations:ATF4_122"
#[25] "deviations:EBF1_67"    "deviations:SREBF1_22"
markerMotifs <- grep("z:", markerMotifs, value = TRUE)#去除掉deviations
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]#去除掉不想分析的
markerMotifs
p <- plotGroups(ArchRProj = projHeme5, 
  groupBy = "Clusters2", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation_single", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)#这是每个图一页，比较直观
pdf("Plot-Groups-Deviations-w-Imputation_single2.pdf",width = 30, height = 10)###这是在一张图，好像没有单独的直观
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
dev.off()
p <- plotEmbedding(
    ArchRProj = projHeme5, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme5)
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
plotPDF(p, name = "Plot-umap_motifmatrix", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)#这是每个图一页，比较直观
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA
################
p <- plotEmbedding(
    ArchRProj = projHeme5, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme5)
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
pdf("Plot-umap-genescorematrix.pdf",width = 20, height = 20)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "Plot-umap-genescorematrix", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)#这是每个图一页，比较直观
#--------------------------FootPrinting----------------------------
motifPositions <- getPositions(projHeme5)
motifPositions
motifs <- c("Tcf7","TBX21","EOMES","PAX8","EBF1","IRF1","IRF8","CEBPA","ATF4","PAX5","JUNB","KLF4","NR4A1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2",force = TRUE)
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)
#####Subtracting the Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias1",
  addDOC = FALSE,
  smoothWindow = 5
)
#####Dividing by the Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
##Co-accessibility with ArchR
projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)
cA
metadata(cA)[[1]]
cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)
cA[[1]]
#cA <- getCoAccessibility(
#    ArchRProj = projHeme5,
#    corCutOff = 0.5,
#    resolution = 1000,#增加resolution会使识别到的GRanges对象中的总条目数减少
#    returnLoops = TRUE
#)
#cA[[1]]
##PLOT——TRACK
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )
#markerGenes <- c("TCF7","TBX21","EOMES","PAX8","EBF1","IRF1","IRF8","CEBPA","ATF4","PAX5","JUNB","KLF4","NR4A1")
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projHeme5)
)
plotPDF(plotList = p, 
    name = "Bcell_Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 5, height = 5)
##Peak2GeneLinkage with ArchR
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = projHeme5,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)
p2g
#markerGenes <- c("TCF7","TBX21","EOMES","PAX8","EBF1","IRF1","IRF8","CEBPA","ATF4","PAX5","JUNB","KLF4","NR4A1")
markerGenes <- c("FOS","FOSL2","FOSB","JUN","JUNB","JUND","BACH2","BACH1","SPI1","SPIC","JDP2","CREB5","NFE2","IRF4")
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projHeme5)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_DE_TF.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 5, height = 5)
###热图
p <- plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2",k=12)##k默认25，可以自己调整聚类数
plotPDF(p, 
    name = "heatmap-Peak2GeneLinks_harmonyk12.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 10, height = 10)
##Identification of Positive TF-Regulators
seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x], useNames = FALSE)
}) %>% Reduce("cbind", .) %>% rowMaxs(useNames = FALSE)
##Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression
corGSM_MM <- correlateMatrices(
    ArchRProj = projHeme5,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGSM_MM
corGIM_MM <- correlateMatrices(
    ArchRProj = projHeme5,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGIM_MM
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
#Step 4. Identify Positive TF Regulators
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
# [1] "BCL11A"    "CEBPA-DT"  "CEBPB"     "CEBPD"     "EBF1"      "EOMES"
# [7] "ETS1"      "ETS2"      "FOS"       "FOSL2"     "FOXO4"     "GABPA"
#[13] "JDP2"      "KLF13"     "LEF1"      "MEF2A"     "MEF2C"     "MITF"
#[19] "NFE2"      "NFIA"      "NFYA"      "PAX5"      "POU2F2"    "PPARD"
#[25] "RARA"      "RUNX1-IT1" "RUNX3"     "SIX5"      "SMAD1"     "SP2"
#[31] "SP4"       "SPI1"      "SPIB"      "STAT6"     "TCF7"      "TP73"
#[37] "ZEB1"
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p
plotPDF(p, 
    name = "GSM_positiveTF.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 5, height = 5)
###用Geneintegrationmatrix替代genescorematrix进行同样运算
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])##没有GSM数量多
##[1] "BACH1"  "CEBPA"  "CEBPB"  "CEBPD"  "EOMES"  "FOS"    "FOSL1"  "ID3"
# [9] "JDP2"   "LEF1"   "MEF2A"  "MEF2C"  "NFE2"   "NFE2L2" "RUNX3"  "SPI1"
#[17] "TCF3"
##----------------------------Myeloid Trajectory - Monocyte Differentiation--------------
p1 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
trajectory <- c("BIN","PC")
trajectory
projHeme5 <- addTrajectory(
    ArchRProj = projHeme5, 
    name = "BIN2PC", 
    groupBy = "Clusters2",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)p <- plotTrajectory(projHeme5, trajectory = "BIN2PC", colorBy = "cellColData", name = "BIN2PC")
p[[1]]
plotPDF(p, name = "Plot-BIN2PC-Traj-UMAP.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 5, height = 5)
head(projHeme5$BIN2PC[!is.na(projHeme5$BIN2PC)])
projHeme5 <- addImputeWeights(projHeme5)
##按某个基因表达做轨迹图
#p1 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneScoreMatrix", name = "NR4A1", continuousSet = "horizonExtra")
#p2 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneIntegrationMatrix", name = "NR4A1", continuousSet = "blueYellow")
#plotPDF(p1, name = "Plot-MyeloidU-Traj-UMAP_NR4A1.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 5, height = 5)
####轨迹Heatmap
trajMM  <- getTrajectory(ArchRProj = projHeme5, name = "BIN2PC", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = projHeme5, name = "BIN2PC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = projHeme5, name = "BIN2PC", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = projHeme5, name = "BIN2PC", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "Plot-BIN2PC-Traj-Heatmaps.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 6, height = 8)
###Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
p <- ht1 + ht2
plotPDF(p, name = "Plot-B2ABC-Traj-Heatmaps_motif.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 10, height = 8)
###
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
p <- ht1 + ht2
plotPDF(p, name = "Plot-B2ABC-Traj-Heatmaps_motif_GIM.pdf", ArchRProj = projHeme5, addDOC = FALSE, width = 10, height = 8)
