

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
scriptPath <- "~/snATAC/B/ArchR/NG_hair_code_ArchR/code/scScalpChromatin-main"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)
# set working directory (The directory of the full preprocessed archr project)
#wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_ProjHeme5_2"
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_MonoSub_ProjHeme5" ###最新分组
plotDir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/p2gLink_plots"
# Color Maps
scriptPath_color <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
allcolour_atac <- readRDS(paste0(scriptPath_color, "/allcolour_atac.rds")) %>% unlist()
allcolour_RNA <- readRDS(paste0(scriptPath_color, "/allcolour_RNA.rds")) %>% unlist()
broadClustCmap <- readRDS(paste0(scriptPath_color, "/broadClustCmap.rds")) %>% unlist()
atacOrder <- c(
 # T/NK / T-cells
  "CD4 NC",
  "CD4 ET", 
  "CD8 NC", 
  "CD8 ET", 
  "NK",
  "NKR", 
  # B/Plasma
  "B IN", # 
  "B Mem",
  "ABC",
  "Plasma",
  # Myeloid
  "Mono C", 
  "Mono NC-I", 
  "Mono NC",
  "Neu",
  "DC"
)
###分亚组
subClusterGroups <- list(
  "T/NK" = c("CD4 NC","CD4 ET","CD8 NC","CD8 ET","NK","NKR"), 
  "Myeloid" = c("Mono C","Mono NC","Mono NC-I","Neu","DC"),
  "Bcells" = c("B IN","B Mem","ABC","Plasma")
  )
subclustered_projects <- c("T/NK", "Myeloid", "Bcells")
# Prepare full-project peak to gene linkages, loops, and coaccessibility (full and subproject links)
##########################################################################################
# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["pbmc"]] <- getPeak2GeneLinks(atac_proj, corCutOff=corrCutoff, resolution = 100)[[1]]
coaccessibility_list <- list()
coAccPeaks <- getCoAccessibility(atac_proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
coAccPeaks$source <- "pbmc"
coaccessibility_list[["pbmc"]] <- coAccPeaks
peak2gene_list <- list()
p2gGR <- getP2G_GR(atac_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
p2gGR$source <- "pbmc"
peak2gene_list[["pbmc"]] <- p2gGR
subclustered_projects <- c("T/NK", "Myeloid","Bcells")
# Retrieve information from subclustered objects
for(subgroup in subclustered_projects){
  message(sprintf("Reading in subcluster %s", subgroup))
  # Read in subclustered project
  sub_dir <- sprintf("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2G_%s", subgroup)
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
idxATAC <- peak2gene_list[["pbmc"]]$idxATAC
names(idxATAC) <- peak2gene_list[["pbmc"]]$peakName
full_p2gGR$idxATAC <- idxATAC[full_p2gGR$peakName]
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing"
# Save lists of p2g objects, etc.
saveRDS(full_p2gGR, file=paste0(wd, "/monoSub_allpeak_multilevel_p2gGR.rds")) # NOT merged or correlation filtered
saveRDS(full_coaccessibility, file=paste0(wd, "/monoSub_allpeak_multilevel_coaccessibility.rds"))
saveRDS(plot_loop_list, file=paste0(wd, "/monoSub_allpeak_multilevel_plot_loops.rds"))
##这里是将整体peak2gene和每个subgroup的peak2gene合起来的总体p2gGR
##########################################################################################
# Upset plot of number of peak to gene links identified per group
##########################################################################################
library(UpSetR)

groups <- unique(full_p2gGR$source)
sub_p2gGR <- full_p2gGR[!is.na(full_p2gGR$Correlation)]

upset_list <- lapply(groups, function(g){
  gr <- sub_p2gGR[sub_p2gGR$source == g & 
    sub_p2gGR$Correlation > corrCutoff & 
    sub_p2gGR$VarQATAC > varCutoffATAC & 
    sub_p2gGR$VarQRNA > varCutoffRNA]
  paste0(gr$peakName, "-", gr$symbol)
  })
names(upset_list) <- groups


plotUpset <- function(plist, main.bar.color="red", keep.order=FALSE){
  # Function to plot upset plots
  upset(
    fromList(plist), 
    sets=names(plist), 
    order.by="freq",
    empty.intersections="off",
    point.size=3, line.size=1.5, matrix.dot.alpha=0.5,
    main.bar.color=main.bar.color,
    keep.order=keep.order
    )
}

pdf(paste0(plotDir, "/Monosub_allpeak_upset_plot_linked_peaks.pdf"), width=10, height=5)
plotUpset(upset_list, main.bar.color="royalblue1", keep.order=FALSE)
dev.off()

##########################################################################################
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
# Plot some comparisons between linked peaks and unlinked peaks
##########################################################################################
p2gGR <- getP2G_GR(atac_proj, corrCutoff=corrCutoff)
#85063
all_linked_peaks <- p2gGR$peakName %>% unique()
#52457
unlinked_peaks <- allPeaksGR[allPeaksGR$peakName %ni% all_linked_peaks]$peakName
#190180
#Most peaks (190180, 78.4%) were not linked to any gene, consistent with the expected small effect size of most CREs32.
df <- data.frame(
  peakName=c(all_linked_peaks, unlinked_peaks),
  label=c(rep("linked", length(all_linked_peaks)), rep("unlinked", length(unlinked_peaks)))
)

df$GC <- allPeaksGR[df$peakName]$GC
df$log10GeneDist <- log10(allPeaksGR[df$peakName]$distToGeneStart)

# Violin plot of GC content between types
wtest <- wilcox.test(df$GC[df$label == "linked"], df$GC[df$label == "unlinked"])
lmedian <- median(df$GC[df$label == "linked"])
ulmedian <- median(df$GC[df$label == "unlinked"])
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)
cmap <- getColorMap(cmaps_BOR$sambaNight, n=2, type="quantitative")
p <- (
  ggplot(df, aes(x=label, y=GC, fill=label))
  + geom_violin(aes(fill=label), adjust = 1.0, scale='width', position=dodge)
  + stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..),
    width=0.75, position=dodge, show.legend = FALSE)
  + scale_color_manual(values=cmap)
  + scale_fill_manual(values=cmap)
  + guides(fill=guide_legend(title=""), 
    colour=guide_legend(override.aes = list(size=5)))
  + ggtitle(sprintf("Wilcoxon test pval: %s\n linked median: %s, unlinked median: %s", wtest$p.value, lmedian, ulmedian))
  + xlab("")
  + ylab("GC content")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1.67,
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1)) 
)
pdf(paste0(plotDir, "/linked_vs_unlinked_GC_content_MonoSub.pdf"), width=6, height=6)
p
dev.off()

# ECDF plot of distance to nearest TSS
p <- (
  ggplot(df, aes(x=log10GeneDist, color=label))
  + stat_ecdf(size=2)
  + geom_vline(xintercept=log10(250000), linetype="dashed") # Dashed line indicating longest possible p2g link
  + scale_color_manual(values=cmap)
  + xlab("Distance to TSS")
  + ylab("Fraction of peaks")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 1)) 
)
pdf(paste0(plotDir, "/linked_vs_unlinked_distance_to_TSS_MonoSub.pdf"), width=6, height=5)
p
dev.off()


##########################################################################################
# Stacked bar chart of types of peaks in p2g links
##########################################################################################
all_peak_types <- getFreqs(allPeaksGR$peakType)
linked_peak_types <- getFreqs(allPeaksGR[allPeaksGR$peakName %in% p2gGR$peakName]$peakType)
multi_peaks <- rownames(g2p_rank_df[g2p_rank_df$ngenes > 1,])
multi_peak_types <- getFreqs(allPeaksGR[allPeaksGR$peakName %in% multi_peaks]$peakType)

plot_df <- data.frame(all_peak_types, linked_peak_types, multi_peak_types)
plot_df <- apply(plot_df, 2, function(x) x/sum(x))
melt_df <- reshape2::melt(plot_df)
colnames(melt_df) <- c("peak_type", "peak_group", "value")

pdf(paste0(plotDir, "/peak_types_barplot_MonoSub.pdf"), width=5, height=6)
stackedBarPlot(melt_df, xcol=2, fillcol=1, ycol=3, 
  cmap=getColorMap(cmaps_BOR$sambaNight, n=4, type="quantitative"), barwidth=0.9)
dev.off()

##########################################################################################
##########################################################################################

# Identify 'highly regulated' genes and 'highly-regulating' peaks
##########################################################################################
# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25
# Coaccessibility cutoffs
coAccCorrCutoff <- 0.4  # Default in getCoAccessibility is 0.5
# Get all peaks
allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
load("~/scRNA/B/10HC11LN/3score_labled15.8immune.combinedCelltypeGroup.RData")
# Get all expressed genes:
count.mat <- Seurat::GetAssayData(object=immune.combined, slot="counts")
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
#p2gFreqs <- data.frame(npeaks=p2gFreqs)##后续有需要关注gene的linkedPeak的number
#write.table(p2gFreqs,file="~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Peak2Gene_rank_df.txt",sep="\t",quote=F,col.name=TRUE)
# Cutoff for defining highly regulated genes (HRGs) determined using elbow rule
cutoff <- 20

# Save HRGs as table
hrg_df <- rank_df[rank_df$npeaks > cutoff,]
hrg_df$gene <- rownames(hrg_df)
table_dir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/supplemental_tables"
write.table(hrg_df, file=paste0(table_dir, "/MonoSub_HRG_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 

# Plot barplot of how many linked peaks per gene
thresh <- 30
threshNpeaks <- rank_df$npeaks
threshNpeaks[threshNpeaks>thresh] <- thresh
nLinkedPeaks <- getFreqs(threshNpeaks)

df <- data.frame(nLinkedPeaks=as.integer(names(nLinkedPeaks)), nGenes=nLinkedPeaks)

pdf(paste0(plotDir, "/MonoSub_nGenes_with_nLinkedPeaks.pdf"), width=8, height=6)
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

pdf(paste0(plotDir, "/Monosub_nPeaks_with_nLinkedGenes.pdf"), width=5, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(g2p_rank_df$npeaks), linetype="dashed")
dev.off()

##########################################################################################
##########################################################################################
# Plot Peak2Gene heatmap
##########################################################################################
################################
for (nclust in c(10,15,20)){
#nclust <- 15
p <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corrCutoff, 
  groupBy="FineClust", 
  nPlot = 1000000, returnMatrices=FALSE, 
  k=nclust, seed=1, palGroup=allcolour_atac
  )
pdf(paste0(plotDir, sprintf("/peakToGeneHeatmap_LabelFineClust_k%s_MonoSub0.45.pdf",nclust)), width=16, height=10)
print(p)
dev.off()
################################
# Need to force it to plot all peaks if you want to match the labeling when you 'returnMatrices'.
p2gMat <- plotPeak2GeneHeatmap(
  atac_proj, 
  corCutOff = corrCutoff, 
  groupBy="FineClust",
  nPlot = 1000000, returnMatrices=TRUE, 
  k=nclust, seed=1)

# Get association of peaks to clusters
kclust_df <- data.frame(
  kclust=p2gMat$ATAC$kmeansId,
  peakName=p2gMat$Peak2GeneLinks$peak,
  gene=p2gMat$Peak2GeneLinks$gene
  )

# Fix peakname
kclust_df$peakName <- sapply(kclust_df$peakName, function(x) strsplit(x, ":|-")[[1]] %>% paste(.,collapse="_"))

# Get motif matches
matches <- getMatches(atac_proj, "Motif")
r1 <- SummarizedExperiment::rowRanges(matches)
rownames(matches) <- paste(seqnames(r1),start(r1),end(r1),sep="_")
matches <- matches[names(allPeaksGR)]

clusters <- unique(kclust_df$kclust) %>% sort()

enrichList <- lapply(clusters, function(x){
  cPeaks <- kclust_df[kclust_df$kclust == x,]$peakName %>% unique()
  ArchR:::.computeEnrichment(matches, which(names(allPeaksGR) %in% cPeaks), seq_len(nrow(matches)))
  }) %>% SimpleList
names(enrichList) <- clusters

# Format output to match ArchR's enrichment output
assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
        enrichList[[y]][colnames(matches),x,drop=FALSE]
      }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
names(assays) <- colnames(enrichList[[1]])
assays <- rev(assays)
res <- SummarizedExperiment::SummarizedExperiment(assays=assays)

formatEnrichMat <- function(mat, topN, minSig, clustCols=TRUE){
  plotFactors <- lapply(colnames(mat), function(x){
    ord <- mat[order(mat[,x], decreasing=TRUE),]
    ord <- ord[ord[,x]>minSig,]
    rownames(head(ord, n=topN))
  }) %>% do.call(c,.) %>% unique()
  pMat <- mat[plotFactors,]
  prettyOrderMat(pMat, clusterCols=clustCols)$mat
}

pMat <- formatEnrichMat(assays(res)$mlog10Padj, 5, 10, clustCols=FALSE)
# Save maximum enrichment
tfs <- strsplit(rownames(pMat), "_") %>% sapply(., `[`, 1)
rownames(pMat) <- paste0(tfs, " (", apply(pMat, 1, function(x) floor(max(x))), ")")

pMat <- apply(pMat, 1, function(x) x/max(x)) %>% t()

pdf(paste0(plotDir, sprintf("/enrichedMotifs_k%s_p2gHM_MonoSub0.45.pdf",nclust)), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pMat, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  #top_annotation = ta,
  row_names_side = "left",
  width = ncol(pMat)*unit(0.5, "cm"),
  height = nrow(pMat)*unit(0.33, "cm"),
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]"
  )
draw(hm)
dev.off()

# GO enrichments of top N genes per cluster 
# ("Top" genes are defined as having the most peak-to-gene links)
source(paste0(scriptPath, "/GO_wrappers.R"))

kclust <- unique(kclust_df$kclust) %>% sort()
all_genes <- kclust_df$gene %>% unique() %>% sort()

# Save table of top linked genes per kclust
nGOgenes <- 250 ##200
topKclustGenes <- lapply(kclust, function(k){
  kclust_df[kclust_df$kclust == k,]$gene %>% getFreqs() %>% head(nGOgenes) %>% names()
  }) %>% do.call(cbind,.)
outfile <- paste0(plotDir, sprintf("/MonoSub0.45_finaltop250_genes_kclust_k%s.tsv", nclust))
write.table(topKclustGenes, file=outfile, quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

GOresults <- lapply(kclust, function(k){
  message(sprintf("Running GO enrichments on k cluster %s...", k))
  clust_genes <- topKclustGenes[,k]
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP") 
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
    )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  })

names(GOresults) <- paste0("cluster_", kclust)

# Plots of GO term enrichments:
pdf(paste0(plotDir, sprintf("/Monosub_kclust_GO_3termsBPonlyBarLim_k%s_250gene0.45.pdf", nclust)), width=10, height=2.5)
for(name in names(GOresults)){
    goRes <- GOresults[[name]]
    if(nrow(goRes)>1){
      print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
        nterms=3, border_color="black", 
        barwidth=0.85, title=name, barLimits=c(0, 15)))
    }
}
dev.off()
}
#################################
# Rank plot:
rank_df$color <- ifelse(rank_df$npeaks >= cutoff, "royalblue1", "black")

# Label select super enhancer - associated genes
label_genes <- c(
"JDP2","ZEB2","CEBPB","FOS","BCL11B","BACH2","TCF7",
"KLF10","RUNX3","ETS1", "IFI30",
"KLF4","IRF8","TCL1A","BLK","PAX5","SPI1","EOMES","PRDM1"
)#peak linked topgenes前20

rank_df$label <- ifelse(rownames(rank_df) %in% label_genes, rownames(rank_df), "")
p <- (ggplot(data=rank_df, aes(x=rank, y=npeaks))
  + ggrepel::geom_text_repel(
          data = rank_df[rank_df$label != "",], aes(x=rank, y=npeaks, label=label), 
          size = 3,
          nudge_x = 2,
          #direction = "x",
          hjust = "outward",
          segment.size = 0.1,
          box.padding=0.5,
          min.segment.length = 0, # draw all segments
          max.overlaps = Inf, # draw all labels
          color = "black")
  + geom_point_rast(aes(color=color))
  + scale_color_manual(values=c("black", "royalblue1"))
  + theme_BOR(border=FALSE)
  + ylab("N linked peaks")
  + xlab("")
#  + ggtitle(sprintf("SE gene enrichment in top %s genes \n -log10(p-value) = %s", k, round(SEmLog10pval,2)))
  + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 1.0,
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
)

pdf(paste0(plotDir, "/nLinkedPeaksPerGene_rastr_Monosub.pdf"))
p
dev.off()
# Assess if 'highly regulated genes' are enriched for super enhancer linked genes
###################################################################################################
# Hnisz super enhancers 2013
young_se_files <- list.files(
  path="~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts/superEnhancers/Adams2015Nature_SEs",
  pattern="*.csv$",
  full.names = TRUE
  )
cell_source <- str_replace(basename(young_se_files), "\\.csv$", "")
names(young_se_files) <- cell_source

young_se_dt <- lapply(names(young_se_files), function(x){
  dt <- fread(young_se_files[x], header=FALSE, sep=",", skip=1)
  dt$source <- x
  dt
  }) %>% rbindlist()
colnames(young_se_dt) <- c("enhID", "chr", "start", "end", "refseq", "enhRank", "SupEnh", "H3K27acDens", "rpmperbp", "source")
young_se_dt <- young_se_dt[SupEnh == 1] # Restrict to only super enhancers

# Convert refseq IDs to gene symbols
library(org.Hs.eg.db)
library(biomaRt)

mart <- useMart("ensembl","hsapiens_gene_ensembl")###需要联网，在登陆节点进行
refseq_to_symbol <- biomaRt::select(
  org.Hs.eg.db, 
  keys=unique(young_se_dt$refseq), 
  columns=c("REFSEQ", "SYMBOL"), 
  keytype="REFSEQ"
  )
ref_to_sym <- refseq_to_symbol$SYMBOL
names(ref_to_sym) <- refseq_to_symbol$REFSEQ
young_se_dt$symbol <- ref_to_sym[young_se_dt$refseq]
young_se_dt <- young_se_dt[!is.na(young_se_dt$symbol)]


#top_SEs <- c(top_SEs_young, top_SEs_fuchs)
top_SEs <- top_SEs_young

all_top_SEs <- unlist(top_SEs) %>% unique()

# Label genes as being SE-associated or not
rank_df$SE <- ifelse(rownames(rank_df) %in% all_top_SEs, 1, 0)

# Hypergeometric enrichment of SE-associated genes in highly-regulated genes
q <- sum(rank_df[rank_df$npeaks > cutoff,]$SE)  # q = number of white balls drawn without replacement
m <- length(all_top_SEs)                        # m = number of white balls in urn
n <- nrow(rank_df) - m                          # n = number of black balls in urn
k <- nrow(rank_df[rank_df$npeaks > cutoff,])    # k = number of balls drawn from urn
SEmLog10pval <- -phyper(q, m, n, k, lower.tail=FALSE, log.p=TRUE)/log(10)

# Get p-value for all sources:
all_sources <- names(top_SEs)
pvals <- sapply(all_sources, function(s){
  rank_df$SE <- ifelse(rownames(rank_df) %in% top_SEs[[s]], 1, 0)
  # Hypergeometric enrichment of SE-associated genes in highly-regulated genes
  q <- sum(rank_df[rank_df$npeaks > cutoff,]$SE)  
  m <- length(top_SEs[[s]])                       
  n <- nrow(rank_df) - m                          
  k <- nrow(rank_df[rank_df$npeaks > cutoff,])    
  -phyper(q, m, n, k, lower.tail=FALSE, log.p=TRUE)/log(10)
  }) %>% sort() %>% rev()

# Plot enrichment of SE's from different sources:
plot_df <- data.frame(rank=1:length(pvals), mlog10padj=-log10(p.adjust(10**-pvals, method="fdr")))
rownames(plot_df) <- names(pvals)

# Save table
write.table(plot_df, file=paste0(table_dir, "/young_SE_enrichment_pval_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 

to_label <- c("CD20","BI_CD4p_CD25-_Il17-_PMAstim_Th", "BI_CD4_Naive_Primary_8pool", "CD19_primary","CD14", "CD3",
 "UCSD_Lung", "BI_Brain_Hippocampus_Middle", "HeLa", "BI_Adipose_Nuclei","UCSD_Spleen", "CD34_fetal", "HepG2")
plot_df$label <- ifelse(rownames(plot_df) %in% to_label, rownames(plot_df), "")
plot_df$color <- "royalblue1"

p <- (ggplot(data=plot_df, aes(x=rank, y=mlog10padj))
  + ggrepel::geom_text_repel(
          data = plot_df[plot_df$label != "",], aes(x=rank, y=mlog10padj, label=label), 
          size = 3,
          nudge_x = 2,
          hjust = "outward",
          segment.size = 0.1,
          box.padding=0.5,
          min.segment.length = 0, # draw all segments
          max.overlaps = Inf, # draw all labels
          color = "black")
  + geom_point(aes(color=color), size=2)
  + scale_color_manual(values=c("royalblue1"))
  + theme_BOR(border=FALSE)
  + ylab("Hypergeometric Enrichment -log10(Padj)")
  + xlab("")
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1.0,
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
)

pdf(paste0(plotDir, "/youngSE_enrichment_by_source.pdf"))
p
dev.off()

# Rank plot:
rank_df$color <- ifelse(rank_df$npeaks >= cutoff, "royalblue1", "black")

# Label select super enhancer - associated genes
label_genes <- c(
  "ICOS", "RUNX3", "TWIST2", "CD84", "CTLA4", "KRT14", "IKZF1", "COL1A1",
  "CD28", "EGFR", "CD3D", "CD69", "ITGAX", "CXCR6",
  "TNF", "RUNX1", "MITF", "FOSL2", "FZD7", "POU2F3", "NR4A1", "CD34"
  )

rank_df$label <- ifelse(rownames(rank_df) %in% label_genes, rownames(rank_df), "")

p <- (ggplot(data=rank_df, aes(x=rank, y=npeaks))
  + ggrepel::geom_text_repel(
          data = rank_df[rank_df$label != "",], aes(x=rank, y=npeaks, label=label), 
          size = 3,
          nudge_x = 2,
          #direction = "x",
          hjust = "outward",
          segment.size = 0.1,
          box.padding=0.5,
          min.segment.length = 0, # draw all segments
          max.overlaps = Inf, # draw all labels
          color = "black")
  + geom_point_rast(aes(color=color))
  + scale_color_manual(values=c("black", "royalblue1"))
  + theme_BOR(border=FALSE)
  + ylab("N linked peaks")
  + xlab("")
  + ggtitle(sprintf("SE gene enrichment in top %s genes \n -log10(p-value) = %s", k, round(SEmLog10pval,2)))
  + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 1.0,
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
)

pdf(paste0(plotDir, "/nLinkedPeaksPerGene_rastr.pdf"))
p
dev.off()
