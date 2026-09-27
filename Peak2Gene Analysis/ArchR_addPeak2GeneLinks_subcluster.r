
###########以Bcells为例######
library(ArchR)
library(dplyr)
library(tidyr)
library(ggrastr)
library(BSgenome.Hsapiens.UCSC.hg19)
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset"
setwd(wd)
full_dir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_MonoSub_ProjHeme5"
subgroups <- "Bcells"
subclustered_projects <- "Bcells"
for(subgroup in subclustered_projects){
  message(sprintf("Reading in subcluster %s", subgroup))
  # Read in subclustered project
  sub_dir <- sprintf("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_ProjHeme5_%s", subgroup)
  sub_proj <- loadArchRProject(sub_dir, force=TRUE)
# Compute group coverages
projHeme4 <- addGroupCoverages(ArchRProj = sub_proj, minCells = 50,
groupBy = "FineClust",force=TRUE)
# The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
# Get peaks that were called on this subproject's subclusters from full ArchR project
full_proj <- loadArchRProject(full_dir, force=TRUE)
full_peaks <- getPeakSet(full_proj)
peaks <- getClusterPeaks(full_proj, clusterNames=unique(projHeme4$FineClust), peakGR=full_peaks)
rm(full_proj)
##
# Now add these peaks to the subproject and generate peak matrix
projHeme4 <- addPeakSet(projHeme4, peakSet=peaks, force=TRUE)
projHeme4 <- addPeakMatrix(projHeme4, force=TRUE)
saveArchRProject(projHeme4)
projHeme4 <- addMotifAnnotations(projHeme4, motifSet="cisbp", name="Motif", force=TRUE)
###
getAvailableMatrices(projHeme4)
#Motif Enrichment in Differential Peaks
plotVarDev <- getVarDeviations(projHeme4, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = paste0(wd,sprintf("Variable-Motif-Deviation-Scores_%s.pdf",subgroup)), width = 5, height = 5, ArchRProj = projHeme4, addDOC = FALSE)
#Co-accessibility with ArchR
projHeme4 <- addCoAccessibility(
    ArchRProj = projHeme4,
    reducedDims = "Harmony"
)
projHeme4 <- addBgdPeaks(projHeme4,force = TRUE)
#Peak2GeneLinkage with ArchR
projHeme4 <- addPeak2GeneLinks(
    ArchRProj = projHeme4,
    reducedDims = "Harmony"
)
#Motif Deviations
projHeme4 <- addDeviationsMatrix(
  ArchRProj = projHeme4, 
  peakAnnotation = "Motif",
  force = TRUE
)
saveArchRProject(ArchRProj = projHeme4, outputDirectory = paste0(wd,sprintf("/peak2Glinks_%s", subgroup)), load = FALSE)
}