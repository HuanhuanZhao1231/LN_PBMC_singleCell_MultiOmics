#!/usr/bin/env Rscript
# 03_LODO_P2G.R
# For each held-out donor:
#   1) PBMC-level P2G on the 8-donor project
#   2) lineage-focused P2G in Lymphoid/Myeloid/Bcells (subGroup的addPeak2GeneLinks之前的dim参数设为15，现在改为默认的30重新计算)
#   3) merge all four sources and deduplicate peak-gene pairs
#   4) apply the same final P2G extraction threshold (correlation >= 0.5)

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(dplyr)
})
#Load Genome Annotations
data("geneAnnoHg19")
data("genomeAnnoHg19")
geneAnno <- geneAnnoHg19
genomeAnno <- genomeAnnoHg19

scriptPath <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source("/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/LODO_P2G/LODO_P2G_scripts/addPeak2Gene_function.R")

# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.4  

cfg <- readRDS(
  "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G/LODO_config.rds"
)

for (heldout in cfg$DONORS) {

  message("\n========================================")
  message("P2G for donor-out: ", heldout)
  message("========================================")

  parent_dir <- file.path(cfg$BASE_OUT, paste0("donor_out_", heldout), "preprocessed")
  proj <- loadArchRProject(parent_dir, force = TRUE)

  # ------------------------------------------------------------
  # A. PBMC-level P2G
  # ------------------------------------------------------------
  message("Running PBMC-level P2G...")

  proj <- addBgdPeaks(proj, force = TRUE)
  metadata(proj@peakSet)$Peak2GeneLinks <- NULL

  proj <- addPeak2GeneLinks_P(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix",
    corCutOff = 0.75,
    dimsToUse = cfg$DIMS_P2G,
    k = cfg$P2G_K,
    knnIteration = cfg$P2G_KNN_ITERATION,
    overlapCutoff = cfg$P2G_OVERLAP_CUTOFF,
    maxDist = cfg$P2G_MAX_DIST,
    scaleTo = cfg$P2G_SCALE_TO,
    log2Norm = cfg$P2G_LOG2_NORM,
    predictionCutoff = cfg$P2G_PREDICTION_CUTOFF,
    addEmpiricalPval = FALSE,
    seed = 1,
    threads = 1,
    verbose = TRUE
  )

  pbmc_p2g <- getP2G_GR(
    proj,
    corrCutoff = NULL,
    varCutoffATAC = -Inf,
    varCutoffRNA = -Inf,
    filtNA = FALSE
  )
  pbmc_p2g$source <- "pbmc"

  # ------------------------------------------------------------
  # B. Lineage-focused P2G
  # ------------------------------------------------------------
  p2g_list <- list(pbmc = pbmc_p2g)

  for (sg in names(cfg$SUBGROUPS)) {

    message("Running subgroup P2G: ", sg)

#    subgroup_cells <- getCellNames(proj)[
#      as.character(proj@cellColData[["FineClust"]]) %in% cfg$SUBGROUPS[[sg]]
#    ]

#    if (length(subgroup_cells) < 100) {
#      warning("Too few cells for subgroup ", sg, " in ", heldout,
#              "; skipping subgroup.")
#      next
#    }

    sg_dir <- file.path(
      cfg$BASE_OUT, paste0("donor_out_", heldout), paste0("subproject_", sg)
    )
sg_proj <- loadArchRProject(sg_dir, force = TRUE)
#     Recreate a subgroup project from the donor-excluded parent.
#    sg_proj <- subsetArchRProject(
#      ArchRProj = proj,
#      cells = subgroup_cells,
#      outputDirectory = sg_dir,
#      dropCells = TRUE,
#      force = TRUE
#    )

    # Recreate group coverages and the fixed, full-cohort-derived subgroup peak set.
    sg_proj <- addGroupCoverages(
      ArchRProj = sg_proj,
      minCells = 50,
      groupBy = "FineClust",
      force = TRUE
    )

    sg_proj <- addPeakSet(
      ArchRProj = sg_proj,
      peakSet = cfg$SUBGROUP_PEAKS[[sg]],
      force = TRUE
    )

    sg_proj <- addPeakMatrix(
      ArchRProj = sg_proj,
      force = TRUE
    )

    sg_proj <- addMotifAnnotations(
      ArchRProj = sg_proj,
      motifSet = "cisbp",
      name = "Motif",
      force = TRUE
    )

    sg_proj <- addCoAccessibility(
      ArchRProj = sg_proj,
      reducedDims = "IterativeLSI"
    )

    sg_proj <- addBgdPeaks(
      ArchRProj = sg_proj,
      force = TRUE
    )
    metadata(sg_proj@peakSet)$Peak2GeneLinks <- NULL

    sg_proj <- addPeak2GeneLinks_P(
      ArchRProj = sg_proj,
      reducedDims = "Harmony"
      #dimsToUse = cfg$DIMS_P2G###上个版本，这里用了full proj的参数，这里需要更正
    )

    sg_p2g <- getP2G_GR(
      sg_proj,
      corrCutoff = NULL,
      varCutoffATAC = -Inf,
      varCutoffRNA = -Inf,
      filtNA = FALSE
    )
    sg_p2g$source <- sg

    p2g_list[[sg]] <- sg_p2g

    saveRDS(
      sg_p2g,
      file.path(cfg$BASE_OUT, paste0("donor_out_", heldout),
                paste0("raw_p2g_addP_", sg, ".rds"))
    )

#    saveArchRProject(
#      ArchRProj = sg_proj,
#      outputDirectory = sg_dir,
#      load = FALSE
#    )

    rm(sg_proj, sg_p2g)
    gc()
  }

  # ------------------------------------------------------------
  # C. Merge + deduplicate
  # ------------------------------------------------------------
  full_p2g <- as(p2g_list, "GRangesList") |> unlist()

  full_p2g <- full_p2g[order(full_p2g$Correlation, decreasing = TRUE)]

  pair_id <- paste0(full_p2g$peakName, "_", full_p2g$symbol)
  filt_p2g <- full_p2g[!duplicated(pair_id)]

  # Keep the same core P2G metadata columns as the original workflow.
  # This intentionally retains the highest-correlation record for each pair.
  p2g_df <- mcols(filt_p2g)[, 1:7]

  metadata(proj@peakSet)$Peak2GeneLinks <- p2g_df

  final_p2g <- getP2G_GR(
    proj,
    corrCutoff = cfg$COR_CUTOFF
  )

  outdir <- file.path(cfg$BASE_OUT, paste0("donor_out_", heldout))
  saveRDS(
    final_p2g,
    file.path(outdir, "final_P2G_addP.rds")
  )
  saveRDS(
    full_p2g,
    file.path(outdir, "merged_raw_P2G_addP.rds")
  )

  # Save a compact summary.
  summary_df <- data.frame(
    heldout = heldout,
    n_cells = nCells(proj),
    n_raw_merged = length(full_p2g),
    n_final = length(final_p2g),
    n_peaks_final = length(unique(final_p2g$peakName)),
    n_genes_final = length(unique(final_p2g$symbol))
  )
  write.csv(
    summary_df,
    file.path(outdir, "P2G_summary_addP.csv"),
    row.names = FALSE
  )

  # Save parent project after P2G calculation.
#  saveArchRProject(
#    ArchRProj = proj,
#    outputDirectory = parent_dir,
#    load = FALSE
#  )

  rm(proj, p2g_list, full_p2g, filt_p2g, final_p2g)
  gc()
}

message("\nAll LODO P2G analyses completed.")
