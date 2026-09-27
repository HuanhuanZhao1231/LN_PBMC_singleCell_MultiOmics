########## Links Fine-mapped SNPs to candidate genes using Peak-to-Gene links
###Links GWAS Significant SNPs to candidate genes与此类似，P < 5e-5
#!/usr/bin/env Rscript

##########################################################################
# Analysis using finemapped SNPs
##########################################################################

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

# Get additional functions, etc.:
scriptPath <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# Set Threads to be used
addArchRThreads(threads = 8)

# set working directory (The directory of the full preprocessed archr project)
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/fine_clustered"
fm_dir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/resources/PICS2"
#Set/Create Working Directory to Folder

##########################################################################################
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
full_p2gGR <- readRDS(file=paste0(p2gGR.dir, "/p2G_allpeak_multilevel_p2gGR.rds")) # NOT merged or correlation filtered
full_coaccessibility <- readRDS(file=paste0(p2gGR.dir, "/p2G_allpeak_multilevel_coaccessibility.rds"))
plot_loop_list <- readRDS(file=paste0(p2gGR.dir, "/p2G_allpeak_multilevel_plot_loops.rds"))

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

# Heatmaps of linked genes and expression by cell type
source(paste0(scriptPath, "/cluster_labels.R"))
#rna_proj$LFineClust <- unlist(rna.FineClust)[rna_proj$FineClust]
#exclude_clust <- c("Unknown", "Cyc.Tc", "Plasma_contam", "TCR.macs", "McSC")
#rna_proj <- rna_proj[,rna_proj$LFineClust %ni% exclude_clust]

countMat <- GetAssayData(object=rna_proj, slot="counts")
groupedCountMat <- averageExpr(countMat, rna_proj$celltype_ABC) # Calculates average log2cp10k, so must only subset afterwards

# Specify order of clusters (Fine Clust)
rnaOrder <- c(
 # Lymphoid / T-cells
 #"Prolif",
  "CD4NC",
  "CD4ET", 
  "CD8NC", 
  "CD8ET", 
  "NK",
  "NKR", 
  # B/Plasma
  "BIN", # 
  "BMem",
  "ABC",
  "Plasma",
  # Myeloid
  "MonoC", 
  "MonoNCI", 
  "MonoNC",
  "cDC",
  "Neu"
  #"LDG",
  #,"pDC",
  #"Mega"
)
rnaOrder <- rnaOrder[rnaOrder %in% colnames(groupedCountMat)]

# Get colors for cluster annotation
# (Link FineClusts to BroadClust cmap)
#cM <- as.matrix(confusionMatrix(rna_proj$LFineClust, rna_proj$BroadClust))
#map_colors <- apply(cM, 1, function(x)colnames(cM)[which.max(x)])
#plotColors <- map_colors[unlist(rna.FineClust)[rnaOrder]]
###统计plotGene的linked_peak的totalNumber需要
rank_df <- read.table("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Peak2Gene_rank_df.txt",header=TRUE)
plotFMsnpsToGenes <- function(expandFmGR, groupedCountMat, trait, cmap, plotDir, backgroundGenes, clusterOrder=NULL, topN=80){
  # Plot finemapped SNPs to genes through peak to gene linkages
  #############################################################
  # expandFmGR = 'expanded' GR of finemapped SNPs linked to genes
  # groupedCountMat = matrix of expression grouped by cell types
  # trait = the GWAS trait/disease to investigate
  # cmap = color map of clusters in groupedCountMat
  # plotDir = plot directory path
  # backgroundGenes = list of genes to use as background for GO term enrichment
  # clusterOrder = order of clusters for heatmap (if null, will bicluster heatmap)
  # topN = number of top genes (by cumulative finemapping probability)
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
  # Cumulative fine-mapping posterior probability per gene
  totalFMprobs <- trait_gr %>% as.data.frame() %>% group_by(linkedGene) %>% summarize(total_FM_prob=sum(fm_prob)) %>% as.data.frame()
  namedFMprobs <- totalFMprobs$total_FM_prob
  names(namedFMprobs) <- totalFMprobs$linkedGene
  message(sprintf("Trait %s has a total of %s unique fmGWAS-linked genes", trait, length(namedFMprobs)))
  # Keep only top N genes for plotting:
  plotGenes <- totalFMprobs[order(totalFMprobs$total_FM_prob, decreasing=TRUE),]$linkedGene %>% head(topN)
  # GO term enrichments
  traitGO <- rbind(
    calcTopGo(backgroundGenes, 
      interestingGenes=plotGenes,
      ontology="MF"),
    calcTopGo(backgroundGenes, 
      interestingGenes=plotGenes,
      ontology="BP")
    )
  traitGO <- traitGO[order(as.numeric(traitGO$pvalue), decreasing=FALSE),]
  pdf(paste0(plotDir, sprintf("/%s_GO_term_enrichments_MFBP_top%s_p2G.pdf", trait_name, topN)), width=8, height=6)
  print(topGObarPlot(traitGO, cmap = cmaps_BOR$comet, 
          nterms=6, border_color="black", 
          barwidth=0.85, title=sprintf("%s GWAS linked Genes", trait)))
  dev.off()
  avgMat <- groupedCountMat[plotGenes,] 
  avgMat <- t(scale(t(avgMat)))
  # Cluster for heatmap
  if(is.null(clusterOrder)){
    plotMat <- prettyOrderMat(avgMat[plotGenes,], clusterCols=TRUE, cutOff=1)$mat
  }else{
    plotMat <- prettyOrderMat(avgMat[plotGenes,clusterOrder], clusterCols=FALSE, cutOff=1)$mat
  }
  plotMat[plotMat > 3] <- 3
  plotMat[plotMat < -3] <- -3
  #colnames(plotMat) <- unlist(rna.FineClust)[colnames(plotMat)]
  # Barplot for number of linked fine-mapped SNPs per gene 
  linkedSNPFreq <- getFreqs(trait_gr$linkedGene)###SNP-Gene不存在重复，那么一个基因有几行就代表有几个snp
  # Barplot for number of linked peaks per gene
  rank_df <- rank_df[plotGenes,,drop=FALSE]
  linkedPeakFreq <- setNames(rank_df$npeaks, rownames(rank_df))
  pdf(paste0(plotDir, sprintf("/%s_GWASlinkedGenesHeatmap_p2G.pdf", trait_name)), width=15, height=12)
  fontsize <- 6
  ht_opt$simple_anno_size <- unit(0.25, "cm")
  ta <- HeatmapAnnotation(
    rna_cluster=rnaOrder,col=list(rna_cluster=allcolour_RNA), 
    show_legend=c(rna_cluster=FALSE), show_annotation_name=c(rna_cluster=FALSE))
  ra <- HeatmapAnnotation(
    nSNPs=anno_barplot(
      linkedSNPFreq[rownames(plotMat)], 
      bar_width=1, height=unit(3, "cm"),
      gp=gpar(fill="red", fontsize=fontsize),
      labels_gp=gpar(fontsize=fontsize)
    ), 
    nPeaks=anno_barplot(
      linkedPeakFreq[rownames(plotMat)], 
      bar_width=1, height=unit(3, "cm"),
      gp=gpar(fill="grey", fontsize=fontsize),
      labels_gp=gpar(fontsize=fontsize)
    ), 
    FMP=anno_barplot(
      namedFMprobs[rownames(plotMat)],
      bar_width=1, height=unit(3, "cm"),
      gp=gpar(fill="blue", fontsize=fontsize), ylim=c(0, max(1, max(namedFMprobs)))
    ),
    which="row")
  hm <- BORHeatmap(
    plotMat, 
    limits=c(-2.0,2.0), 
    clusterCols=FALSE, clusterRows=FALSE,
    labelCols=TRUE, labelRows=TRUE,
    dataColors = cmaps_BOR$sunrise,
    top_annotation = ta,
    right_annotation = ra,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    width = ncol(plotMat)*unit(0.4, "cm"),
    height = nrow(plotMat)*unit(0.25, "cm"),
    legendTitle="row Z-score",
    border_gp = gpar(col="black") # Add a black border to entire heatmap
    )
  draw(hm)
  dev.off()
}
# Plot selected traits
plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Systemic lupus erythematosus", 
#cmap=plotColors, 
plotDir=plotDir, backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)

plotFMsnpsToGenes(expandFmGR, groupedCountMat, "Estimated glomerular filtration rate",
 #cmap=plotColors, 
 plotDir=plotDir, 
  backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)#28Gene

plotFMsnpsToGenes(expandFmGR, groupedCountMat, "eGFR", 
#cmap=plotColors, 
plotDir=plotDir, backgroundGenes=unique(p2gGR$symbol), clusterOrder=rnaOrder, topN=topN)


############################################################################################################

# Enrichment of TFs for AGA

trait_gr <- expandFmGR[expandFmGR$disease_trait %in% c("Systemic lupus erythematosus")]
trait_gr <- trait_gr[order(trait_gr$fm_prob, decreasing=TRUE)]
# Further filter by fine-mapping posterior probability
trait_gr <- trait_gr[trait_gr$fm_prob >= 0.01]
trait_gr <- trait_gr[!is.na(trait_gr$linkedGene)]
trait_gr <- trait_gr[!duplicated(trait_gr$SNP_to_gene)]

# Cumulative fine-mapping posterior probability per gene
totalFMprobs <- trait_gr %>% as.data.frame() %>% group_by(linkedGene) %>% summarize(total_FM_prob=sum(fm_prob)) %>% as.data.frame()
namedFMprobs <- totalFMprobs$total_FM_prob
names(namedFMprobs) <- totalFMprobs$linkedGene

# Keep only top N genes for plotting:
plotGenes <- totalFMprobs[order(totalFMprobs$total_FM_prob, decreasing=TRUE),]$linkedGene %>% head(topN)

trait_genes <- trait_gr$linkedGene %>% unique()
bg_genes <- unique(p2gGR$symbol)
TFdb <- fread("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/resources/lambert_2018_TF_master_list.txt")
TFs <- TFdb %>% filter(Is_TF == "Yes") %>% pull(Symbol)
valid_TFs <- TFs[TFs %in% bg_genes]

# Hypergeometric enrichment
q <- sum(trait_genes %in% valid_TFs)  # q = number of white balls drawn without replacement
m <- length(valid_TFs)                # m = number of white balls in urn
n <- length(bg_genes) -  m            # n = number of black balls in urn
k <- length(trait_genes)              # k = number of balls drawn from urn
TFphypPval <- phyper(q, m, n, k, lower.tail=FALSE, log.p=FALSE)

enrichment <- (q/length(trait_genes))/(m/length(bg_genes))
OR <- (q/(k-q))/(m/n)

# 查看结果
print(result)
###################################################################################



dotPlot <- function(df, xcol, ycol, color_col, size_col, xorder=NULL, yorder=NULL, cmap=NULL,
  color_label=NULL, size_label=NULL, aspectRatio=NULL, sizeLims=NULL, colorLims=NULL){
  # Plot rectangular dot plot where color and size map to some values in df
  # (Assumes xcol, ycol, color_col and size_col are named columns)

  # If neither x or y col order is provided, make something up
  # Sort df:
  if(is.null(xorder)){
    xorder <- unique(df[,xcol]) %>% sort()
  }
  if(is.null(yorder)){
    yorder <- unique(df[,ycol]) %>% sort()
  }
  if(is.null(aspectRatio)){
    aspectRatio <- length(yorder)/length(xorder) # What is the best aspect ratio for this chart?
  }
  df[,xcol] <- factor(df[,xcol], levels=xorder)
  df[,ycol] <- factor(df[,ycol], levels=yorder)
  df <- df[order(df[,xcol], df[,ycol]),]

  # Make plot:
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], color=df[,color_col], size=ifelse(df[,size_col] > 0, df[,size_col], NA)))
    + geom_point()
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(),
            panel.grid.minor= element_blank(),
            plot.margin = unit(c(0.25,0,0.25,1), "cm"),
            aspect.ratio = aspectRatio,
            axis.text.x = element_text(angle = 90, hjust = 1))
    + theme(axis.text.y = element_blank())###加了这句不显示y轴的标签
    + guides(
      fill = guide_legend(title=""),
      colour = guide_legend(title=color_label, override.aes = list(size=5)),
      size = guide_legend(title=size_label)
      )
  )
  if(!is.null(cmap)){
    if(!is.null(colorLims)){
      p <- p + scale_color_gradientn(colors=cmap, limits=colorLims, oob=scales::squish, name = "")
    }else{
      p <- p + scale_color_gradientn(colors=cmap, name = "")
    }
  }
  if(!is.null(sizeLims)){
    p <- p + scale_size_continuous(limits=sizeLims)
  }
  p
}
