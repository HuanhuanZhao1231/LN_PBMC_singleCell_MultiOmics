library(FigR)
library(ArchR)
library(igraph)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggrastr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(circlize)
library(networkD3)
library(network)
library(tibble)
HRG <- read.table("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/supplemental_tables/rmUn_p2G_HRG_table.tsv",header=TRUE,sep="\t")
dorcGenes <- HRG$gene
dorcMat <- getDORCScores(ATAC.se = peakmatrix, # Has to be same SE as used in previous step
                         dorcTab = mycisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 1)
########################smooth these (sparse) DORC counts#####################
setwd("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/FigR/data")
# Get cell KNNs
lsi_mat <- getReducedDims(
    ArchRProj = proj,
    reducedDims = "IterativeLSI")

cellkNN <- FNN::get.knn(lsi_mat,k=30)$nn.index
rownames(cellkNN) <- colnames(dorcMat)
# Smooth dorc scores using cell KNNs (k=30)
library(doParallel)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN,mat = dorcMat,nCores = 8)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
colnames(RNAmat_matched) <- colnames(peakmatrix) # Just so that the smoothing function doesn't throw an error (matching cell barcodes in the KNN and the matrix)
RNAmat.s <- smoothScoresNN(NNmat = cellkNN,mat = RNAmat_matched,nCores = 8)
save(dorcMat.s,RNAmat.s,file="SmoothMat_dorcMat.s_RNAmat.s.RData")
#########################TF-gene associations############################
#As before, we can now determine TF-gene associations, and begin inferring a regulatory network based on the previous DORC definitions, together with information drawn from a database of TF binding sequence motifs
figR.d <- runFigRGRN(ATAC.se = peakmatrix, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = mycisCorr.filt, # Filtered peak-gene associations
                     genome = "hg19",
                     dorcMat = dorcMat.s,
                     dorcK = 5, 
                     rnaMat = RNAmat.s, 
                     nCores = 1)
save(figR.d,file="allHRG_figR.d.RData")
library(ComplexHeatmap)
# 添加 Status 标记
figR.d$Status <- "notSig"  # 默认标记为 notSig
figR.d$Status[figR.d$DORC %in% HRG_pos$X] <- "up"
figR.d$Status[figR.d$DORC %in% HRG_neg$X] <- "down"
pdf("allHRG_plotfigRHeatmap_score2_up_down.pdf",height=50,width=5)
plotfigRHeatmap_colour(figR.d = figR.d,
                score.cut = 2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()
pdf("upHRG_plotfigRHeatmap_score2.pdf",height=40,width=5)
plotfigRHeatmap(figR.d = figR.d
                score.cut = 1.5,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()
pdf("downHRG_plotfigRHeatmap_Score2.pdf",height=6.5,width=4)#
plotfigRHeatmap(figR.d = downHRG,
                score.cut = 2,
                cluster_columns = TRUE,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
#Using absolute score cut-off of: 1.5 ..
#Plotting 69 DORCs x 20TFs
dev.off()
#################Network view#######
library(networkD3)
label <- read.table("node_label_downHRG.txt",header=F,sep="\t")#
node_label_downHRG <- label$V1
plotfigRNetwork.2.part(figR.d,
                      score.cut = 2,
                      weight.edges = TRUE,
                      legend=FALSE,##TRUE显示图例
                      charge = -15,
                      fontSize= 12,
                     size_dorc = 10,
                     size_TF = 18,
                     genes.to.label = node_label_downHRG)