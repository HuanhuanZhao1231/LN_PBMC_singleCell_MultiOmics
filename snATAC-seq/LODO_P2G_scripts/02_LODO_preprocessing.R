#!/usr/bin/env Rscript
# 02_LODO_preprocessing.R
# Leave-one-donor-out preprocessing:
# exclude one LN donor -> recompute IterativeLSI -> recompute Harmony.
# The full 242,637-peak universe is retained/fixed for downstream comparability.

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
})
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
CONFIG_RDS <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G/LODO_config.rds"
cfg <- readRDS(CONFIG_RDS)

for (heldout in cfg$DONORS) {
  message("\n==============================")
  message("LODO held-out donor: ", heldout)
  message("==============================")

  outdir <- file.path(cfg$BASE_OUT, paste0("donor_out_", heldout), "preprocessed")
  dir.create(dirname(outdir), recursive = TRUE, showWarnings = FALSE)

  # Do not rebuild if the finished project is already present.
  if (dir.exists(outdir) && file.exists(file.path(outdir, "ArrowFiles"))) {
    message("Project directory exists; loading and checking.")
  }

  full_proj <- loadArchRProject(cfg$FULL_DIR, force = TRUE)

  keep_cells <- cfg$CELL_META$cell[cfg$CELL_META$donor != heldout]
  stopifnot(length(keep_cells) > 0)

  proj <- subsetArchRProject(
    ArchRProj = full_proj,
    cells = keep_cells,
    outputDirectory = outdir,
    dropCells = TRUE,
    force = TRUE
  )

  # Recompute the representation after donor exclusion.
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 3,
    sampleCellsFinal = 50000,
    projectCellsPre = TRUE,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 15000,
    dimsToUse = cfg$DIMS_ATAC,
    force = TRUE
  )

  # Recompute Harmony without the held-out donor.
  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )

  # Sanity checks.
  lsi_dim <- dim(getReducedDims(proj, reducedDims = "IterativeLSI"))
  har_dim <- dim(getReducedDims(proj, reducedDims = "Harmony"))

  if (lsi_dim[2] != length(cfg$DIMS_ATAC)) {
    stop("Unexpected IterativeLSI dimension in ", heldout)
  }
  if (har_dim[2] != length(cfg$DIMS_ATAC)) {
    stop("Unexpected Harmony dimension in ", heldout)
  }

  saveRDS(
    list(
      heldout = heldout,
      n_cells = nCells(proj),
      donors_remaining = sort(unique(as.character(proj@cellColData[["Sample"]]))),
      IterativeLSI_dim = lsi_dim,
      Harmony_dim = har_dim
    ),
    file.path(cfg$BASE_OUT, paste0("donor_out_", heldout), "preprocessing_qc.rds")
  )

  saveArchRProject(
    ArchRProj = proj,
    outputDirectory = outdir,
    load = FALSE
  )

  rm(full_proj, proj)
  gc()
}
