#!/usr/bin/env Rscript
# 04_compare_full_LODO.R
# Compare each LODO network with the full 9-donor network at:
# link level, target-gene level, correlation concordance, and HRG recovery.

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
})
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

cfg <- readRDS(
  "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G/LODO_config.rds"
)

# ------------------------------------------------------------
# 1. Full reference P2G network
# ------------------------------------------------------------
full_proj <- loadArchRProject(cfg$FULL_DIR, force = TRUE)
full_p2gGR <- readRDS(file="~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/p2G_allpeak_multilevel_p2gGR.rds") # NOT merged or correlation filtered
# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(full_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)
# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))] %>% sort()
# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(full_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF
# Get full merged p2g links
full_p2g <- getP2G_GR(full_proj, corrCutoff = cfg$COR_CUTOFF)



make_pair_id <- function(x) {
  paste0(x$peakName, "_", x$symbol)
}

full_id <- make_pair_id(full_p2g)
full_genes <- unique(as.character(full_p2g$symbol))

# HRG definition: >=20 P2G links.
HRG_THRESHOLD <- 20
full_gene_counts <- table(as.character(full_p2g$symbol))
full_hrg <- names(full_gene_counts[full_gene_counts > HRG_THRESHOLD])

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
}

results <- lapply(cfg$DONORS, function(heldout) {

  f <- file.path(cfg$BASE_OUT, paste0("donor_out_", heldout), "final_P2G.rds")
  lodo <- readRDS(f)

  lodo_id <- make_pair_id(lodo)
  shared_id <- intersect(full_id, lodo_id)

  # Correlation concordance on shared peak-gene pairs.
  fidx <- match(shared_id, full_id)
  lidx <- match(shared_id, lodo_id)

  full_cor <- as.numeric(full_p2g$Correlation[fidx])
  lodo_cor <- as.numeric(lodo$Correlation[lidx])

  shared_gene <- intersect(
    unique(as.character(full_p2g$symbol)),
    unique(as.character(lodo$symbol))
  )

  lodo_gene_counts <- table(as.character(lodo$symbol))
  lodo_hrg <- names(lodo_gene_counts[lodo_gene_counts > HRG_THRESHOLD])

  data.frame(
    heldout = heldout,
    n_cells = nrow(cfg$CELL_META[cfg$CELL_META$donor != heldout, , drop = FALSE]),
    full_links = length(full_p2g),
    lodo_links = length(lodo),
    shared_links = length(shared_id),
    link_recovery = length(shared_id) / length(full_id),
    link_jaccard = length(shared_id) /
      length(union(full_id, lodo_id)),
    full_genes = length(full_genes),
    lodo_genes = length(unique(as.character(lodo$symbol))),
    shared_genes = length(shared_gene),
    gene_recovery = length(shared_gene) / length(full_genes),
    gene_jaccard = length(shared_gene) /
      length(union(full_genes, unique(as.character(lodo$symbol)))),
    Pearson_r = safe_cor(full_cor, lodo_cor, "pearson"),
    Spearman_rho = safe_cor(full_cor, lodo_cor, "spearman"),
    full_HRG = length(full_hrg),
    lodo_HRG = length(lodo_hrg),
    shared_HRG = length(intersect(full_hrg, lodo_hrg)),
    HRG_recovery = length(intersect(full_hrg, lodo_hrg)) / length(full_hrg),
    HRG_jaccard = length(intersect(full_hrg, lodo_hrg)) /
      length(union(full_hrg, lodo_hrg))
  )
}) |> bind_rows()

write.csv(
  results,
  file.path(cfg$BASE_OUT, "04_full_vs_LODO_summary.csv"),
  row.names = FALSE
)

# Save shared-link correlation tables for downstream plotting.
shared_tables <- lapply(cfg$DONORS, function(heldout) {
  lodo <- readRDS(file.path(
    cfg$BASE_OUT, paste0("donor_out_", heldout), "final_P2G.rds"
  ))
  full_id <- make_pair_id(full_p2g)
  lodo_id <- make_pair_id(lodo)
  shared_id <- intersect(full_id, lodo_id)

  data.frame(
    heldout = heldout,
    pair_id = shared_id,
    full_Correlation = as.numeric(full_p2g$Correlation[match(shared_id, full_id)]),
    LODO_Correlation = as.numeric(lodo$Correlation[match(shared_id, lodo_id)])
  )
}) |> bind_rows()

write.csv(
  shared_tables,
  file.path(cfg$BASE_OUT, "04_shared_link_correlations.csv"),
  row.names = FALSE
)

print(results)

message("\nComparison completed.")
