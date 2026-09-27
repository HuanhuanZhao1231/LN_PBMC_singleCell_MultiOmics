library(ArchR)
library(dplyr)
library(tidyr)
library(ggrastr)
library(BSgenome.Hsapiens.UCSC.hg19)
setwd("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset")
proj_LN9 <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_Monosub_ProjHeme4")
projHeme4 <- addGroupCoverages(ArchRProj = proj_LN9, minCells = 50,
groupBy = "FineClust",force=TRUE)
# The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
pathToMacs2 <- "~/.conda/envs/macs2/bin/macs2"
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "FineClust",
    peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
    pathToMacs2 = pathToMacs2,
    force=TRUE
)
getPeakSet(projHeme4)
saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-proj_N9_Monosub_ProjHeme4", load = FALSE)
projHeme5 <- addPeakMatrix(projHeme4)
getAvailableMatrices(projHeme5)
#Motif Enrichment in Differential Peaks
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
#Motif Deviations
projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
#Co-accessibility with ArchR
projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)
#Peak2GeneLinkage with ArchR
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = projHeme5,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = TRUE
)
saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-proj_N9_Monosub_ProjHeme5", load = FALSE)
