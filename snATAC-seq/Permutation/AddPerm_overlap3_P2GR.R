library(ArchR)
library(igraph)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggrastr)
#Load Genome Annotations
data("geneAnnoHg19")
data("genomeAnnoHg19")
geneAnno <- geneAnnoHg19
genomeAnno <- genomeAnnoHg19
# Get additional functions, etc.:##很重要，下面会用到文件里的自定义function
scriptPath <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/perm_functions.R"))
# Set Threads to be used
#addArchRThreads(threads = 16)
# set working directory (The directory of the full preprocessed archr project)
plotDir <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2gLink_plots"
# Color Maps
scriptPath_color <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
allcolour_atac <- readRDS(paste0(scriptPath_color, "/allcolour_atac.rds")) %>% unlist()
allcolour_RNA <- readRDS(paste0(scriptPath_color, "/allcolour_RNA.rds")) %>% unlist()
broadClustCmap <- readRDS(paste0(scriptPath_color, "/broadClustCmap.rds")) %>% unlist()
##########
#atac_proj <- loadArchRProject("/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_rmUn_MonoSub_ProjHeme5")
#perm_proj <- addPermPeak2GeneLinks(
#    ArchRProj=atac_proj,
#    reducedDims="IterativeLSI",
#    corCutOff=0.5,
#    overlapCutoff = 0.3,
#    addPermutedPval=TRUE,
#    nperm=1000,
#    seed=123,
#    threads=16
#)

#p2g_perm <- metadata(perm_proj@peakSet)$Peak2GeneLinks


#summary(p2g_perm$PermFDR)


#sum(
#p2g_perm$Correlation>0.5 &
##na.rm=TRUE
#)
setwd("/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset")
#saveArchRProject(
#    ArchRProj = perm_proj,
#    outputDirectory="Save-Perm_OC3_proj_N9_rmUn_MonoSub_ProjHeme5_before_P2G",
#    load=FALSE
#)
#########################################重新加载perm_proj########

subclustered_projects <- c(
"Lymphoid",
"Myeloid",
"Bcells"
)

################3个subGroup重新addpermPeak2Genelinks##########
run_subproject_permP2G <- function(
    project_dir,
    name,
    save_dir,
    nperm=1000
){


    message(
        "Loading ",
        name
    )


    proj <- loadArchRProject(
        project_dir,
        force=TRUE
    )


    message(
        "Running permutation P2G for ",
        name
    )


    proj <- addPermPeak2GeneLinks(
        ArchRProj = proj,
        reducedDims = "IterativeLSI",
        corCutOff = 0.5,
        overlapCutoff = 0.3,
        addPermutedPval = TRUE,
        nperm = nperm,
        seed = 123,
        threads = 8
    )


    saveArchRProject(
        ArchRProj = proj,
        outputDirectory = save_dir,
        load = FALSE
    )


    return(
        list(
            proj=proj
        )
    )

}
#lymph_result <- run_subproject_permP2G(
    
#    project_dir =
#    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_Lymphoid",
    
#    name="Lymphoid",
    
#    save_dir =
#    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC3_peak2Glinks_Lymphoid",
    
#    nperm=1000
#)
myeloid_result <- run_subproject_permP2G(
    
    project_dir =
    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_Myeloid",
    
    name="Myeloid",
    
    save_dir =
    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC3_peak2Glinks_Myeloid",
    
    nperm=1000
)
Bcells_result <- run_subproject_permP2G(
    
    project_dir =
    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_Bcells",
    
    name="Bcells",
    
    save_dir =
    "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC3_peak2Glinks_Bcells",
    
    nperm=1000
)

