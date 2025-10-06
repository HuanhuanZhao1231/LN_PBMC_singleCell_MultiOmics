library(ArchR)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggrastr)
# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "~/snATAC/B/ArchR/Save-ProjHeme3_2_11_LSI3dim15"
plotDir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/p2gLink_plots"
setwd(wd)

#Load Genome Annotations
data("geneAnnoHg19")
data("genomeAnnoHg19")
geneAnno <- geneAnnoHg19
genomeAnno <- genomeAnnoHg19

atac_proj <- loadArchRProject(wd, force=TRUE)
rna_proj_path <- "/public/home/zhaohuanhuan/snATAC/B/ArchR/18type_imm.sce.rds"

# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.4  # Default in getCoAccessibility is 0.5

subClusterGroups <- list(
  "T/NK" = c("aTc1","aTc2","aTc3","aTc4","aTc5","aTc6","aTc7"), 
  "Myeloid" = c("aMy1","aMy2","aMy3","aMy4","aMy5","aMy6"),
  "Bcells" = c("aBc1","aBc2")
  )
subClusterCells <- lapply(subClusterGroups, function(x){
  getCellNames(proj)[as.character(proj@cellColData[["BroadClust"]]) %in% x]
  })

subClusterArchR <- function(proj, subCells, outdir){
  # Subset an ArchR project for focused analysis

  message(sprintf("Subgroup has %s cells.", length(subCells)))
  sub_proj <- subsetArchRProject(
      ArchRProj = proj,
      cells = subCells,
      outputDirectory = outdir,
      dropCells = TRUE,
      force = TRUE
  )}
saveArchRProject(sub_proj)

# Get all peaks
allPeaksGR <- getPeakSet(ArchRProj=atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

##########################################################################################
# Prepare full-project peak to gene linkages, loops, and coaccessibility (full and subproject links)
##########################################################################################
subclustered_projects <- c("T/NK", "Myeloid", "Bcells")

# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["PBMC"]] <- getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff, resolution = 100)[[1]]
coaccessibility_list <- list()
coAccPeaks <- getCoAccessibility(atac_proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
coAccPeaks$source <- "PBMC"
coaccessibility_list[["PBMC"]] <- coAccPeaks
peak2gene_list <- list()
p2gGR <- getP2G_GR(atac_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
p2gGR$source <- "PBMC"
peak2gene_list[["PBMC"]] <- p2gGR
# Retrieve information from subclustered objects
for(subgroup in subclustered_projects){
  message(sprintf("Reading in subcluster %s", subgroup))
  # Read in subclustered project
  sub_dir <- sprintf("~/scATAC_preprocessing/subclustered_%s", subgroup)
  sub_proj <- loadArchRProject(sub_dir, force=TRUE)

  # Get sub-project p2g links
  subP2G <- getP2G_GR(sub_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
  subP2G$source <- subgroup
  peak2gene_list[[subgroup]] <- subP2G

  # Get sub-project loops
  plot_loop_list[[subgroup]] <- getPeak2GeneLinks(sub_proj, corCutOff=corrCutoff, resolution = 100)[[1]]

  # Get coaccessibility
  coAccPeaks <- getCoAccessibility(sub_proj, corCutOff=coAccCorrCutoff, returnLoops=TRUE)[[1]]
  coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
  coAccPeaks$source <- subgroup
  coaccessibility_list[[subgroup]] <- coAccPeaks
}

full_p2gGR <- as(peak2gene_list, "GRangesList") %>% unlist()
full_coaccessibility <- as(coaccessibility_list, "GRangesList") %>% unlist()

# Fix idxATAC to match the full peak set
idxATAC <- peak2gene_list[["scalp"]]$idxATAC
names(idxATAC) <- peak2gene_list[["scalp"]]$peakName
full_p2gGR$idxATAC <- idxATAC[full_p2gGR$peakName]

# Save lists of p2g objects, etc.
saveRDS(full_p2gGR, file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
saveRDS(full_coaccessibility, file=paste0(wd, "/multilevel_coaccessibility.rds"))
saveRDS(plot_loop_list, file=paste0(wd, "/multilevel_plot_loops.rds"))
# Filter redundant peak to gene links
##########################################################################################
# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(atac_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(atac_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

# Identify 'highly regulated' genes and 'highly-regulating' peaks
##########################################################################################

# Get all expressed genes:
count.mat <- Seurat::GetAssayData(object=readRDS(rna_proj_path), slot="counts")
minUMIs <- 1
minCells <- 2
valid.genes <- rownames(count.mat[rowSums(count.mat > minUMIs) > minCells,])

# Get distribution of peaks to gene linkages and identify 'highly-regulated' genes
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)
p2gFreqs <- getFreqs(p2gGR$symbol)
valid.genes <- c(valid.genes, unique(p2gGR$symbol)) %>% unique()

noLinks <- valid.genes[valid.genes %ni% names(p2gFreqs)]
zilch <- rep(0, length(noLinks))
names(zilch) <- noLinks
p2gFreqs <- c(p2gFreqs, zilch)
x <- 1:length(p2gFreqs)
rank_df <- data.frame(npeaks=p2gFreqs, rank=x)

# Cutoff for defining highly regulated genes (HRGs) determined using elbow rule
cutoff <- 20

# Save HRGs as table
hrg_df <- rank_df[rank_df$npeaks > cutoff,]
hrg_df$gene <- rownames(hrg_df)
table_dir <- "/oak/stanford/groups/wjg/boberrey/hairATAC/results/supplemental_tables"
write.table(hrg_df, file=paste0(table_dir, "/HRG_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 

# Plot barplot of how many linked peaks per gene
thresh <- 30
threshNpeaks <- rank_df$npeaks
threshNpeaks[threshNpeaks>thresh] <- thresh
nLinkedPeaks <- getFreqs(threshNpeaks)

df <- data.frame(nLinkedPeaks=as.integer(names(nLinkedPeaks)), nGenes=nLinkedPeaks)

pdf(paste0(plotDir, "/nGenes_with_nLinkedPeaks.pdf"), width=8, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(rank_df$npeaks), linetype="dashed")
dev.off()

# Get distribution of gene to peak linkages 
valid.peaks <- allPeaksGR$peakName
g2pFreqs <- getFreqs(p2gGR$peakName)

noLinks <- valid.peaks[valid.peaks %ni% names(g2pFreqs)]
zilch <- rep(0, length(noLinks))
names(zilch) <- noLinks
g2pFreqs <- c(g2pFreqs, zilch)
x <- 1:length(g2pFreqs)
g2p_rank_df <- data.frame(ngenes=g2pFreqs, rank=x)

# Plot barplot of how many linked peaks per gene
thresh <- 5
threshNgenes <- g2p_rank_df$ngenes
threshNgenes[threshNgenes>thresh] <- thresh
nLinkedGenes <- getFreqs(threshNgenes)

df <- data.frame(nLinkedGenes=as.integer(names(nLinkedGenes)), nPeaks=nLinkedGenes)

pdf(paste0(plotDir, "/nPeaks_with_nLinkedGenes.pdf"), width=5, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(g2p_rank_df$npeaks), linetype="dashed")
dev.off()
################################FigR################################
