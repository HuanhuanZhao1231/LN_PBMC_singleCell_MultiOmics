#!/usr/bin/env Rscript
# 05_HRG_reproducibility.R
# Quantify reproducibility of highly regulated genes (HRGs) under LODO.

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cfg <- readRDS(
  "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G/LODO_config.rds"
)
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
HRG_THRESHOLD <- 20

full_proj <- loadArchRProject(cfg$FULL_DIR, force = TRUE)
# Full reference.
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

count_hrg <- function(p2g, threshold = 20) {
  tb <- table(as.character(p2g$symbol))
  names(tb[tb > threshold])
}

full_hrg <- count_hrg(full_p2g, HRG_THRESHOLD)

# HRG membership matrix.
hrg_list <- list(Full = full_hrg)
linkcount_list <- list(
  data.frame(
    gene = names(table(as.character(full_p2g$symbol))),
    Full = as.integer(table(as.character(full_p2g$symbol)))
  )
)

for (heldout in cfg$DONORS) {
  p2g <- readRDS(file.path(
    cfg$BASE_OUT, paste0("donor_out_", heldout), "final_P2G.rds"
  ))

  hrg_list[[paste0("out_", heldout)]] <- count_hrg(p2g, HRG_THRESHOLD)

  tb <- table(as.character(p2g$symbol))
  linkcount_list[[paste0("out_", heldout)]] <- data.frame(
    gene = names(tb),
    count = as.integer(tb)
  )
}

all_genes <- sort(unique(unlist(hrg_list)))
hrg_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

for (nm in names(hrg_list)) {
  hrg_membership[[nm]] <- all_genes %in% hrg_list[[nm]]
}

write.csv(
  hrg_membership,
  file.path(cfg$BASE_OUT, "05_HRG_membership_matrix.csv"),
  row.names = FALSE
)

# Number of LODO iterations in which each gene remains an HRG.
lodo_names <- paste0("out_", cfg$DONORS)
hrg_membership$LODO_n <- rowSums(
  hrg_membership[, lodo_names, drop = FALSE]
)

hrg_membership$Full_HRG <- hrg_membership$Full

write.csv(
  hrg_membership,
  file.path(cfg$BASE_OUT, "05_HRG_reproducibility_by_gene.csv"),
  row.names = FALSE
)

# Recovery summary.
hrg_summary <- bind_rows(lapply(cfg$DONORS, function(heldout) {
  hrg <- hrg_list[[paste0("out_", heldout)]]
  data.frame(
    heldout = heldout,
    full_HRG = length(full_hrg),
    LODO_HRG = length(hrg),
    shared_HRG = length(intersect(full_hrg, hrg)),
    recovery = length(intersect(full_hrg, hrg)) / length(full_hrg),
    jaccard = length(intersect(full_hrg, hrg)) /
      length(union(full_hrg, hrg))
  )
}))

write.csv(
  hrg_summary,
  file.path(cfg$BASE_OUT, "05_HRG_recovery_summary.csv"),
  row.names = FALSE
)

# Stability categories among full HRGs.
full_hrg_repro <- hrg_membership %>%
  filter(Full_HRG) %>%
  arrange(desc(LODO_n), gene)

write.csv(
  full_hrg_repro,
  file.path(cfg$BASE_OUT, "05_full_HRG_stability.csv"),
  row.names = FALSE
)

cat("\nFull HRGs (>= ", HRG_THRESHOLD, " links): ", length(full_hrg), "\n", sep = "")
print(hrg_summary)

cat("\nFull HRGs retained in >=8/9 LODO analyses: ",
    sum(full_hrg_repro$LODO_n >= 8), "\n", sep = "")
cat("Full HRGs retained in all 9 LODO analyses: ",
    sum(full_hrg_repro$LODO_n == 9), "\n", sep = "")

# A simple HRG recovery plot for the supplement.
p <- ggplot(hrg_summary, aes(x = heldout, y = recovery)) +
  geom_col() +
  geom_hline(yintercept = 0.8, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Held-out donor",
    y = "Full-network HRG recovery",
    title = "HRG reproducibility under leave-one-donor-out analysis"
  ) +
  theme_classic()

ggsave(
  file.path(cfg$BASE_OUT, "05_HRG_recovery.pdf"),
  p, width = 6.5, height = 4.5
)
