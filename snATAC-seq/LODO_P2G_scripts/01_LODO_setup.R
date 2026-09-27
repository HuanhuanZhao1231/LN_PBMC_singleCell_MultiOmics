#!/usr/bin/env Rscript
# 01_LODO_setup.R
# Create a reproducible configuration for the 9-donor leave-one-donor-out analysis.
# This script does not modify the original ArchR project.

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(GenomicRanges)
})

set.seed(1)
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

FULL_DIR <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_rmUn_MonoSub_ProjHeme5"
BASE_OUT <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G"
dir.create(BASE_OUT, recursive = TRUE, showWarnings = FALSE)

DONORS <- c("LN1","LN2","LN3","LN4","LN5","LN6","LN9","LN10","LN11")

SUBGROUPS <- list(
  Lymphoid = c("CD4 NC","CD4 ET","CD8 NC","CD8 ET","NK","NKR"),
  Myeloid  = c("Mono C","Mono NC","Mono NC-I","Neu","DC"),
  Bcells   = c("B IN","B Mem","ABC","Plasma")
)

# Current ArchR project has 15 dimensions in both reductions.
DIMS_ATAC <- 1:15
DIMS_P2G <- 1:15

COR_CUTOFF <- 0.5
VAR_CUTOFF_ATAC <- 0.25
VAR_CUTOFF_RNA <- 0.25
COACC_CUTOFF <- 0.4

# Parameters recorded for the original P2G workflow.
P2G_RESOLUTION <- 100
P2G_K <- 100
P2G_KNN_ITERATION <- 500
P2G_OVERLAP_CUTOFF <- 0.8
P2G_MAX_DIST <- 250000
P2G_SCALE_TO <- 10000
P2G_LOG2_NORM <- TRUE
P2G_PREDICTION_CUTOFF <- 0.4

full_proj <- loadArchRProject(FULL_DIR, force = TRUE)

# Confirm the dimensions before any LODO run.
lsi_dim <- dim(getReducedDims(full_proj, reducedDims = "IterativeLSI"))
har_dim <- dim(getReducedDims(full_proj, reducedDims = "Harmony"))
if (lsi_dim[2] != 15 || har_dim[2] != 15) {
  stop("Expected 15 dimensions in both IterativeLSI and Harmony.")
}

CELL_META <- data.frame(
  cell = getCellNames(full_proj),
  donor = as.character(full_proj@cellColData[["Sample"]]),
  FineClust = as.character(full_proj@cellColData[["FineClust"]]),
  stringsAsFactors = FALSE
)
rownames(CELL_META) <- CELL_META$cell

# Fixed full peak set.
FULL_PEAKS <- getPeakSet(full_proj)
FULL_PEAKS$peakName <- paste0(
  as.character(seqnames(FULL_PEAKS)), "_",
  start(FULL_PEAKS), "_", end(FULL_PEAKS)
)
names(FULL_PEAKS) <- FULL_PEAKS$peakName

# Fixed subgroup peak sets derived ONCE from the full 9-donor project.
# This preserves the same peak universe across all LODO iterations while
# retaining the original lineage-focused peak-calling logic.
SUBGROUP_PEAKS <- lapply(names(SUBGROUPS), function(sg) {
  cl <- SUBGROUPS[[sg]]
  present <- intersect(cl, unique(CELL_META$FineClust))
  getClusterPeaks(
    full_proj,
    clusterNames = present,
    peakGR = FULL_PEAKS
  )
})
names(SUBGROUP_PEAKS) <- names(SUBGROUPS)

# Save counts for QC.
donor_counts <- CELL_META %>%
  dplyr::count(donor, name = "n_cells")

subgroup_counts <- lapply(names(SUBGROUPS), function(sg) {
  CELL_META %>%
    dplyr::filter(FineClust %in% SUBGROUPS[[sg]]) %>%
    dplyr::count(FineClust, name = "n_cells") %>%
    dplyr::mutate(Subgroup = sg)
}) %>% dplyr::bind_rows()

cfg <- list(
  FULL_DIR = FULL_DIR,
  BASE_OUT = BASE_OUT,
  DONORS = DONORS,
  SUBGROUPS = SUBGROUPS,
  DIMS_ATAC = DIMS_ATAC,
  DIMS_P2G = DIMS_P2G,
  COR_CUTOFF = COR_CUTOFF,
  VAR_CUTOFF_ATAC = VAR_CUTOFF_ATAC,
  VAR_CUTOFF_RNA = VAR_CUTOFF_RNA,
  COACC_CUTOFF = COACC_CUTOFF,
  P2G_RESOLUTION = P2G_RESOLUTION,
  P2G_K = P2G_K,
  P2G_KNN_ITERATION = P2G_KNN_ITERATION,
  P2G_OVERLAP_CUTOFF = P2G_OVERLAP_CUTOFF,
  P2G_MAX_DIST = P2G_MAX_DIST,
  P2G_SCALE_TO = P2G_SCALE_TO,
  P2G_LOG2_NORM = P2G_LOG2_NORM,
  P2G_PREDICTION_CUTOFF = P2G_PREDICTION_CUTOFF,
  CELL_META = CELL_META,
  FULL_PEAKS = FULL_PEAKS,
  SUBGROUP_PEAKS = SUBGROUP_PEAKS,
  DONOR_COUNTS = donor_counts,
  SUBGROUP_COUNTS = subgroup_counts,
  FULL_LSI_DIM = lsi_dim,
  FULL_HARMONY_DIM = har_dim
)

saveRDS(cfg, file.path(BASE_OUT, "LODO_config.rds"))
write.csv(donor_counts, file.path(BASE_OUT, "full_donor_cell_counts.csv"), row.names = FALSE)
write.csv(subgroup_counts, file.path(BASE_OUT, "full_subgroup_cell_counts.csv"), row.names = FALSE)

cat("\nLODO setup completed.\n")
cat("Full cells:", nrow(CELL_META), "\n")
cat("Full peaks:", length(FULL_PEAKS), "\n")
cat("Config:", file.path(BASE_OUT, "LODO_config.rds"), "\n")
