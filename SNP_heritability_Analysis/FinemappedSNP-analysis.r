#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(plyranges)
  library(data.table)
  library(stringr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(parallel)
  library(ggrepel)
  library(ComplexHeatmap)
})

# Set Threads to be used
ncores <- 8
addArchRThreads(threads = ncores)
# set working directory (The directory of the full preprocessed archr project)
wd <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/fine_clustered"
fm_dir <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/resources/PICS2"
#Set/Create Working Directory to Folder
#dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

topN <- 80 # Number of genes to plot in heatmap

##########################################################################################
# Preparing Data
##########################################################################################
wd <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_rmUn_MonoSub_ProjHeme5"
atac_proj <- loadArchRProject(wd, force=TRUE)
load("/public/home/zhaohuanhuan/scRNA/B/10HC11LN/3score_labled15.8immune.combinedCelltypeGroup_celltype_ABC.RData")
rna_proj <- immune.combined
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")
allcolour_atac <- readRDS(paste0(scriptPath, "/allcolour_atac_noBlankforLDSCplot.rds")) %>% unlist()
allGenesGR <- getGenes(atac_proj)
raw_finemapped_gr <- readRDS(paste0(fm_dir, "/unfiltered_finemapping_genomic_range.rds"))
# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
raw_finemapped_gr <- raw_finemapped_gr[order(raw_finemapped_gr$fm_prob, decreasing=TRUE)]
raw_finemapped_gr$trait_snp <- paste0(raw_finemapped_gr$disease_trait, "_", raw_finemapped_gr$linked_SNP)
raw_finemapped_gr <- raw_finemapped_gr[!duplicated(raw_finemapped_gr$trait_snp)] %>% sort()
# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
###########
disease_traits <- list(
lupus_other_autoimmune=c(
    "Asthma",
    "Cutaneous lupus erythematosus",
    "Systemic lupus erythematosus",
    "Crohn's disease",
    "Rheumatoid arthritis"
  ),
lupus_related_traits=c(
    "Interferon alpha levels in systemic lupus erythematosus",
    "C-reactive protein levels",
    "Complement C3 and C4 levels",
    "Response to cyclophosphamide in systemic lupus erythematosus with lupus nephritis"
),
immune_related_nephropathy=c(
    "IgA nephropathy",
    "Membranous nephropathy"
),
other_kidney_disease=c(
    "Diabetic kidney disease",
    "Chronic kidney disease"
  ),
kidney_related_traits=c(
    "eGFR",##UKBB
    "Urea",#UKBB
    "Urinary albumin-to-creatinine ratio"
),
brain_traits=c(
    "Parkinson's disease",
    "Schizophrenia",
    "Neuroticism"
  ),
  other_traits=c(
    "BMI", # Finacune finemapped##UKBB
    "Height",##UKBB
    "SBP" # Finacune finemapped##UKBB
  ))
all_disease_traits <- do.call(c,disease_traits) %>% unname()

pList <- lapply(all_disease_traits, function(dt){
  message(sprintf("Testing trait %s...", dt))
  ol <- getPeakOLpct(raw_finemapped_gr[raw_finemapped_gr$disease_trait == dt], peaks_gr=allPeaksGR, breaks=prob_bins)
  pct_ol <- ol$ol_pct
  nSNPs <- ol$ol_nSNPs
  xlabs <- paste0(names(nSNPs), sprintf("\n(%s)", nSNPs))
  df <- data.frame(fm_probs=names(pct_ol), percent_snps_overlapping=pct_ol)
  df$fm_probs <- factor(df$fm_probs, levels=names(pct_ol), ordered=TRUE)
  qcBarPlot(df, cmap="grey80", barwidth=0.9, border_color="black") + 
    scale_y_continuous(limits=c(0, 0.4), expand = c(0, 0)) + 
    ggtitle(dt) + scale_x_discrete(labels=xlabs)
  })
pdf("pctPeaksOL_by_fm_probs_eachTrait2.pdf", width=5, height=5)
pList
dev.off()
# Links Fine-mapped SNPs to candidate genes using Peak-to-Gene links
##########################################################################################
finemapped_GR <- readRDS(paste0(fm_dir, "/filtered_finemapping_genomic_range.rds"))
# Some of the fine-mapped SNPs are duplicated (i.e. the Finacune SNPs sometimes have both FINEMAP and SuSiE finemapping results)
# Deduplicate trait-SNP pairs prior to proceeding with enrichment analyses:
finemapped_GR <- finemapped_GR[order(finemapped_GR$fm_prob, decreasing=TRUE)]
finemapped_GR$trait_snp <- paste0(finemapped_GR$disease_trait, "_", finemapped_GR$linked_SNP)
finemapped_GR <- finemapped_GR[!duplicated(finemapped_GR$trait_snp)] %>% sort()

# Load full project p2g links, plot loops, etc.
p2gGR.dir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing"
full_p2gGR <- readRDS(file=paste0(p2gGR.dir, "/rmUn_p2G_allpeak_multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_coaccessibility <- readRDS(file=paste0(p2gGR.dir, "/rmUn_p2G_allpeak_multilevel_coaccessibility.rds"))
plot_loop_list <- readRDS(file=paste0(p2gGR.dir, "/rmUn_p2G_allpeak_multilevel_plot_loops.rds"))

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$symbol, "-", full_p2gGR$peakName))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

##########################################################################################
# Identify Finemapped SNPs linked to genes
##########################################################################################

p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff, 
  varCutoffATAC=varCutoffATAC, varCutoffRNA=varCutoffRNA, filtNA=TRUE)

pol <- findOverlaps(p2gGR, finemapped_GR, type="any", maxgap=-1L, ignore.strand=TRUE)
expandFmGR <- finemapped_GR[to(pol)]
expandFmGR$linkedGene <- p2gGR[from(pol)]$symbol
expandFmGR$linkedPeak <- p2gGR[from(pol)]$peakName
expandFmGR$p2gCorr <- p2gGR[from(pol)]$Correlation
expandFmGR$SNP_to_gene <- paste(expandFmGR$linked_SNP, expandFmGR$linkedGene, sep="_")

trait_name <- strsplit(trait, split=" ")[[1]] %>% paste(., collapse="_")
  trait_gr <- expandFmGR[expandFmGR$disease_trait %in% c(trait)]
  trait_gr <- trait_gr[order(trait_gr$fm_prob, decreasing=TRUE)]
  # Further filter by finemapping posterior probability
  trait_gr <- trait_gr[trait_gr$fm_prob >= 0.01]
  trait_gr <- trait_gr[!is.na(trait_gr$linkedGene)]
  trait_gr <- trait_gr[!duplicated(trait_gr$SNP_to_gene)]
  # Save fmGWAS-linked genes
  message(sprintf("Saving fmGWAS-genes for traint %s...", trait))
  saveRDS(trait_gr, paste0(plotDir, sprintf("/%s_fmGWAS_genes_p2G.rds", trait_name)))