#!/usr/bin/env Rscript
# 06_make_Supplementary_Figure.R
# Make a compact donor-aware P2G reproducibility figure.
# Panels:
# A. Link recovery
# B. Gene recovery
# C. Correlation concordance
# D. HRG recovery
# E. HRG stability across LODO iterations

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))


BASE_OUT <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/LODO_P2G"

res <- read.csv(
  file.path(BASE_OUT, "04_full_vs_LODO_summary.csv"),
  stringsAsFactors = FALSE
)

# Consistent donor ordering.
donor_order <- c("LN1","LN2","LN3","LN4","LN5","LN6","LN9","LN10","LN11")
res$heldout <- factor(res$heldout, levels = donor_order)

# ------------------------------------------------------------
# A: Link recovery
# ------------------------------------------------------------
pA <- ggplot(res, aes(x = heldout, y = link_recovery)) +
  geom_col() +
  geom_hline(yintercept = 0.8, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "A",
    x = "Held-out donor",
    y = "P2G link recovery"
  ) +
  theme_classic(base_size = 11)

# ------------------------------------------------------------
# B: Gene recovery
# ------------------------------------------------------------
pB <- ggplot(res, aes(x = heldout, y = gene_recovery)) +
  geom_col() +
  geom_hline(yintercept = 0.8, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "B",
    x = "Held-out donor",
    y = "Target-gene recovery"
  ) +
  theme_classic(base_size = 11)

# ------------------------------------------------------------
# C: Correlation concordance
# ------------------------------------------------------------
cor_long <- res %>%
  dplyr::select(heldout, Pearson_r, Spearman_rho) %>%
  pivot_longer(
    cols = c(Pearson_r, Spearman_rho),
    names_to = "metric",
    values_to = "value"
  )

pC <- ggplot(cor_long, aes(x = heldout, y = value, group = metric)) +
  geom_point() +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "C",
    x = "Held-out donor",
    y = "Correlation with full network",
    caption = "Pearson and Spearman correlations were calculated for shared peak-gene links."
  ) +
  theme_classic(base_size = 11)

# ------------------------------------------------------------
# D: HRG recovery
# ------------------------------------------------------------
pD <- ggplot(res, aes(x = heldout, y = HRG_recovery)) +
  geom_col() +
  geom_hline(yintercept = 0.8, linetype = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "D",
    x = "Held-out donor",
    y = "HRG recovery"
  ) +
  theme_classic(base_size = 11)

# ------------------------------------------------------------
# E: HRG stability
# ------------------------------------------------------------
hrg <- read.csv(
  file.path(BASE_OUT, "05_full_HRG_stability.csv"),
  stringsAsFactors = FALSE
)

hrg_plot <- hrg %>%
  count(LODO_n, name = "n_HRG") %>%
  complete(LODO_n = 0:9, fill = list(n_HRG = 0))

pE <- ggplot(hrg_plot, aes(x = LODO_n, y = n_HRG)) +
  geom_col() +
  scale_x_continuous(breaks = 0:9) +
  labs(
    title = "E",
    x = "Number of LODO analyses retaining HRG status",
    y = "Number of full-network HRGs"
  ) +
  theme_classic(base_size = 11)

fig <- (pA | pB) / (pC | pD)

ggsave(
  file.path(BASE_OUT, "06_Supplementary_Figure_LODO_P2G.pdf"),
  fig,
  width = 9,
  height = 7
)

ggsave(
  file.path(BASE_OUT, "06_Supplementary_Figure_LODO_P2G.png"),
  fig,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(BASE_OUT, "06_Supplementary_Figure_HRG_stability.pdf"),
  pE,
  width = 6,
  height = 4.5
)

message("Supplementary figures written to: ", BASE_OUT)
