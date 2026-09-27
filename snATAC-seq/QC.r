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