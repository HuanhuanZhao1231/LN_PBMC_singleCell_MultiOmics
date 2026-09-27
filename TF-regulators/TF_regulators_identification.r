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
#library(GGally)
library(network)
library(tibble)
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
# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.4  # Default in getCoAccessibility is 0.5
####
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_MonoSub_ProjHeme5"
proj <- loadArchRProject(wd)

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
wd <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing"
full_p2gGR <- readRDS(file=paste0(wd, "/p2G_allpeak_multilevel_p2gGR.rds"))
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(proj@peakSet)$Peak2GeneLinks <- new_p2g_DF
# Plot some comparisons between linked peaks and unlinked peaks
##########################################################################################
p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff)####p2gGR替代cisCor.filt
#mycisCorr <- as.data.frame(p2gGR)
mycisCorr.filt <- as.data.frame(p2gGR)
mycisCorr.filt <- mycisCorr.filt[,c("idxATAC","peakName","symbol","Correlation","FDR")]
head(mycisCorr.filt)
  idxATAC             peakName   symbol Correlation          FDR
1      40   chr1_935283_935783     HES4   0.5247590 1.950283e-35
2     114 chr1_1143313_1143813 TNFRSF18   0.6941894 2.819511e-71
3     116 chr1_1148087_1148587 TNFRSF18   0.6207895 5.096884e-53
4     122 chr1_1151738_1152238 TNFRSF18   0.5356248 3.793607e-37
5     123 chr1_1152366_1152866 TNFRSF18   0.6387602 4.978996e-57
6     110 chr1_1136628_1137128  TNFRSF4   0.5008834 6.862846e-32
#####将peakName改成chr1:762060-762560格式
mycisCorr.filt$peakName <- gsub("_(\\d+)_(\\d+)", ":\\1-\\2", mycisCorr.filt$peakName)
####列名改成与FigR计算的cisCorr一致
colnames(mycisCorr.filt) <- c("Peak","PeakRanges","Gene","rObs","pvalZ")
save(mycisCorr.filt, file="mycisCorr.filt_fromP2gGR_full_p2gGR_archR.RData")
######主要看HRG,读取HRGlist##########
HRG <- read.table("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/supplemental_tables/p2G_HRG_table.tsv",header=TRUE,sep="\t")
dorcGenes <- HRG$gene
######peak的矩阵####
dorcMat <- getDORCScores(ATAC.se = peakmatrix, # Has to be same SE as used in previous step
                         dorcTab = mycisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 1)
####提取LSI矩阵##
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
######################Ranking TF drivers#####################
rankDrivers(figR.d,rankBy = "meanScore")
#> Warning: ggrepel: 63 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps
ggsave(file="allHRG_rank_drivers_meanScore2.pdf",height=5,width=5)
rankDrivers(figR.d,score.cut = 1.5,rankBy = "nTargets",interactive = FALSE)
ggsave(file="allHRG_rank_drivers_nTargets.pdf",height=5,width=5)
#You can also set the interactive parameter (which the rankDrivers function takes) to TRUE, which will let you hover over points and provide some extra information in terms of the number of estimated activated vs repressed genes
#######################Heatmap view###############################
library(ComplexHeatmap)
# 添加 Status 标记
figR.d$Status <- "notSig"  # 默认标记为 notSig
figR.d$Status[figR.d$DORC %in% HRG_pos$X] <- "up"
figR.d$Status[figR.d$DORC %in% HRG_neg$X] <- "down"
####plotfigRHeatmap_colour自定义function，增加字体颜色，见此页末端######
pdf("allHRG_plotfigRHeatmap_score2_up_down.pdf",height=50,width=5)
plotfigRHeatmap_colour(figR.d = figR.d,
                score.cut = 2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()
######给ISG添加Status标记
######ISGlist##########
ISG <- read.table("~/scRNA/B/10HC11LN/ISG/IRG100genes_27celltypes_SHEHC_fromCell_useISGgenelist.txt",header=TRUE)
ISGgene <- ISG$Gene
allcell_deg <- read.csv("~/scRNA/B/10HC11LN/DEG/data/LNvsHC_EdgeR_DEG_immune.combined.csv",sep=",") 
    allcell_deg.pos <- subset(allcell_deg, logFC >= 0)
    allcell_deg.neg <- subset(allcell_de2, logFC < 0)
ISG_pos <- allcell_deg.pos[allcell_deg.pos$X %in% ISGgene,]
ISG_neg <- allcell_deg.neg[allcell_deg.neg$X %in% ISGgene,]###0个基因
figR.d$Status <- "notSig"  # 默认标记为 notSig
figR.d$Status[figR.d$DORC %in% ISG_pos$X] <- "up"
pdf("ISG_plotfigRHeatmap_score2_up_down.pdf",height=6,width=3.5)
plotfigRHeatmap_colour(figR.d = figR.d,
                score.cut = 2,
                column_names_gp = gpar(fontsize=8), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()
##################upHRG################
pdf("upHRG_plotfigRHeatmap_score2.pdf",height=40,width=5)
plotfigRHeatmap(figR.d = figR.d
                score.cut = 2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
##################单独downHRG################
pdf("ISG_plotfigRHeatmap2.pdf",height=6.5,width=4)
plotfigRHeatmap(figR.d = figR.d,
                score.cut = 2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
#Using absolute score cut-off of: 1.5 ..
#Plotting 68 DORCs x 19TFs
dev.off()
pdf("Coloumcluster_down_plotfigRHeatmap_Score2.pdf",height=6.5,width=4)
plotfigRHeatmap(figR.d = downHRG,
                score.cut = 2,
                cluster_columns = TRUE,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
#Using absolute score cut-off of: 1.5 ..
#Plotting 69 DORCs x 20TFs
dev.off()
########################缩小heatmap################
color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE) 
color_use <- color_peak
res.h <- plotfigRHeatmap(figR.d = figR.d,
                score.cut = 2,
                #DORCs = genes.to.label[genes.to.label %in% figR.d$DORC],
                #TFs = tf.to.label[tf.to.label %in% figR.d$Motif],
                #column_names_gp = gpar(fontsize=8,fontface = 'bold'), # from ComplexHeatmap
                #row_names_gp = gpar(fontsize=10),
                show_row_dend = FALSE # from ComplexHeatmap
                )
res.h.mat <- res.h@matrix
set.seed(123);row_idx <- row_order(res.h)
set.seed(123);col_idx <- column_order(res.h)
rowid_dorc_list <- rownames(res.h.mat)[row_idx]
colid_tf_list <- colnames(res.h.mat)[col_idx]
regMat <- res.h.mat[row_idx,col_idx]
#saveRDS(regMat,'regMat.score_cutoff_2.0.rds')
######基因太多，个别展示###############
genes.to.label_allHRG <- c(
'FOSB','SPI1','NLRP3','RETN','STAT5A','CSF3R',
'S100A8','KLF4','CEBPA','PRKCD','PTPRJ','MYD88','KLF10','IL1B',
'CCRL2','ETS2','IL17RA','IFNAR1','IL10RB','JDP2','FOS',
'ACSL1','TLR2','CCR2','IRF8',
'FOSL2','ZEB2','IFNG',
'TYK2','CD38','IGLL5',
'HLA-DMB','HLA-DOB','HLA-DQA1',
'CCL4','PDCD1',"TBX21",
'TCF7','RUNX2','BACH2',
'KLF2','ORMDL3','TCL1A','RNASET2','FCRL5',
'FAM167A','BLK','ETS1')###包括一些低表达的淋巴细胞基因
genes.to.label <- c(
'TLR2','BCL6','GPR27','IL1B','CCR2','CCRL2','MEF2A',
'ZEB2','FOSL2','HLA-DMB','FTL','CSF2RB',
'IRF8','STAT5A','SPI1','NLRP3','RETN','STAT5A','CSF3R',
'S100A8','KLF4','CEBPA','PRKCD','PTPRJ','MYD88','KLF10','IL1B',
'CCRL2','ETS2','IL17RA','IFNAR1','IL10RB','JDP2','FOS',
'IFNG','TYK2','CD38','IGLL5','PDCD1',
'HLA-DMB','HLA-DRB1','HLA-DPA1',
'TLR1','IFNGR1','PRAM1','LY86','CD1D','CLEC7A',
'LILRA6','LILRA3','LILRA1','HLA-DMA','HLA-DPA1','HLA-DRA',
"IFI30","IRF5","IRF8","ITGAM","ITGAX","NCF2","SIRPB1","SLC16A5","RELT",###SLE-Genetic variant落入p2GR的基因
"FOXO3","VEGFA","SLC43A2","IRF5","FTH1","KLF4","FTL","STAT6","SP1",
"LRP1","ACTR2","IFIT3","SPI1","MXD1","VDR","CD14","SLC11A1","FOSB",
"VASP","NEAT1","TLR2","NRIP1","CTSB"#######eGFR-Genetic variant落入p2GR的基因
)####全是upHRG

> intersect(unique(trait_gr$linkedGene),upHRG$DORC)
 [1] "SCARB2"    "FOXO3"     "LRRC25"    "ZFHX3"     "MCL1"      "ADAMTSL4"
 [7] "CTSS"      "VEGFA"     "SLC43A2"   "RILP"      "SCARF1"    "IRF5"
[13] "BEST1"     "FTH1"      "SSR1"      "RREB1"     "COTL1"     "CRISPLD2"
[19] "RCOR1"     "MAP3K11"   "AP5B1"     "FRMD4B"    "RIN3"      "FBXL19"
[25] "LINC00482" "TET2"      "GNA15"     "ZBTB7B"    "ADAM15"    "UQCRC1"
[31] "PFKFB4"    "NAMPT"     "KLF4"      "ARID3A"    "SBNO2"     "FTL"
[37] "NRIP1"     "CTSD"      "SMARCD3"   "ACVR1B"    "INSR"      "CTSH"
[43] "LRRC8D"    "ZDHHC7"    "SLC16A6"   "SLC38A10"  "MIR22HG"   "LRP1"
[49] "STAT6"     "STAC3"     "TTYH3"     "PER1"      "SP1"       "ACTR2"
[55] "TPCN2"     "RAB1A"     "IFIT3"     "SPI1"      "PPARG"     "ZCCHC24"
[61] "PPIF"      "DOK3"      "RAB24"     "CUX1"      "MXD1"      "RBPJ"
[67] "VDR"       "CD14"      "CTSB"      "SEMA4A"    "LAMTOR2"   "SMG5"
[73] "HBEGF"     "SRA1"      "SLC11A1"   "GPBAR1"    "MXD3"      "FOSB"
[79] "VASP"      "SH3BP2"    "NEAT1"     "CDC42EP2"  "PRKAR1A"   "KIAA0513"
[85] "TXNDC5"    "FCN1"      "CFD"       "RBM47"     "RNF149"    "TLR2"
[91] "RPS6KA4"

idx <- vector()
for(i in genes.to.label){
    idx <- c(idx,grep(pattern = paste0('^',i,'$') , x = rownames(regMat) ) )
    
}

idx.id <- rownames(regMat)[idx]

ha_right <- rowAnnotation(foo = anno_mark(at = idx,  #HeatmapAnnotation(... which = 'row')
                                    labels = idx.id,
                                    #labels = as.character(ID_map.use$Gene),
                                    labels_gp = gpar(fontsize = 8),###右侧基因名称
                                    link_width = unit(2, "mm"),
                                    extend = unit(10, "mm"),
                                    padding = unit(0.8, "mm") #important
                                    #annotation_width= unit(4, "cm")
                                   )

) 

res.hp = Heatmap(#regMat.sel, name = "TF-target", 
                 regMat, name = "Score", 
         cluster_rows = FALSE, 
         cluster_columns = TRUE,
         show_row_names = FALSE,#FALSE,
         show_column_names = TRUE,#FALSE,
         show_row_dend = FALSE,
         show_column_dend = FALSE,
         use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
         column_names_gp = gpar(fontsize = 8),
         ##col = circlize::colorRamp2(seq(-1.5,1.5,by=3/10), viridis(n = 11,option = "C")),
         col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_use)-1)), color_use),
         #col = color_use,
         #col = circlize::colorRamp2(seq(-1,1,by=2/(length(color_gradient_my)-1)), color_gradient_my),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_tfdev1)-1)), color_tfdev1),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/255), color_peak),
         na_col = 'white',
         #column_km = 3,
         #row_km = 3,
         #heatmap_legend_param = list(color_bar = "continuous"),
         #right_annotation = ha#,heatmap_width=unit(8, "cm"),
         ##clustering_distance_rows  = 'pearson', 
         ##clustering_distance_columns  = 'pearson',
         #column_split = cluster[,1],
         column_gap = unit(.1,'cm'),
         row_gap = unit(.1,'cm'),
         #column_labels = levels(cluster[,1]),
         column_names_side = 'bottom',
         #top_annotation = ha_top,#trajectory bin cluster belonging
         #bottom_annotation = , #trajectory arrow and text, no use decorate_heatmap_body
         ##left_annotation = ha_left, #peak dar cluster belonging
         right_annotation = ha_right #peak annotation gene label
)
pdf(file = 'Coloum-cluster-withSLE_eGFR_Geneticgene_upHRG_logFC0.5_TF-mining-regulation.heatmap.withlable.pdf',width = 5,height = 8.5,useDingbats = FALSE,fonts = NULL)
res.hp
dev.off()
#################Network view#############是在线的图，需要借助浏览器，所以保存figr.d去本地Rstudio操作##
library(networkD3)
##默认的plotfigRNetwork的参数不好看所以自定义参数画图
plotfigRNetwork.2 <- function(figR.d,
                              score.cut=1,
                              DORCs=NULL,
                              TFs=NULL,
                              size_dorc = 8,
                              size_TF = 10,
                              charge = -15,
                              legend = FALSE,
                              fontSize = 8,
                              weight.edges=FALSE){
  # Network view
  
  # Filter
  net.dat <- figR.d %>% filter(abs(Score) >= score.cut)
  
  if(!is.null(DORCs))
    net.dat <- net.dat %>% filter(DORC %in% DORCs)
  
  if(!is.null(TFs))
    net.dat <- net.dat %>% filter(Motif %in% TFs)
  
  net.dat$Motif <- paste0(net.dat$Motif, ".")
  net.dat$DORC <- paste0(net.dat$DORC)
  
  dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", size = size_dorc)
  tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", size = size_TF)
  nodes <- rbind(dorcs, tfs)

  # 添加 fontcolor 列，所有文本颜色设为黑色
  nodes$fontcolor <- "black"

  edges <- as.data.frame(net.dat)

  # Make edges into links (subtract 1 for 0 indexing)
  links <- data.frame(source=unlist(lapply(edges$Motif, function(x) {which(nodes$name==x)-1})), 
                      target=unlist(lapply(edges$DORC, function(x) {which(nodes$name==x)-1})), 
                      corr=edges$Corr,
                      enrichment=edges$Enrichment.P)

  links$Value <- scales::rescale(abs(edges$Score)) * 10 + 1

  # Set of colors you can choose from for TF/DORC nodes
  colors <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Tomato", 
              "Forest Green", "Sky Blue", "Gray", "Steelblue3", "Firebrick2", 
              "Brown", "darkgrey")
  nodeColorMap <- data.frame(color = colors, hex = gplots::col2hex(colors))

  getColors <- function(tfColor, dorcColor = NULL) {
    temp <- c(as.character(nodeColorMap[nodeColorMap$color==tfColor,]$hex),
              as.character(nodeColorMap[nodeColorMap$color==dorcColor,]$hex))
    if (is.null(dorcColor)) {
      temp <- temp[1]
    }
    colors <- paste(temp, collapse = '", "')
    colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
    colorJS
  }

  network_plot <- forceNetwork(
  Links = links, 
  Nodes = nodes,
  Source = "target",
  Target = "source",
  NodeID = "name",
  Group = "group",
  Value = "Value",
  Nodesize = "size",
  radiusCalculation = "Math.sqrt(d.nodesize)*2",
  arrows = FALSE,
  opacityNoHover = 0.6,
  opacity = 1,
  zoom = TRUE,
  bounded = TRUE,
  charge = charge,  
  fontSize = fontSize,
  legend = legend,
  fontFamily = "Helvetica",
  colourScale = getColors(tfColor = "Tomato", dorcColor = "darkgrey"),
  linkColour = ifelse(links$corr > 0, as.character(nodeColorMap[nodeColorMap$color=="Brown",]$hex),
                      as.character(nodeColorMap[nodeColorMap$color=="Blue",]$hex)),
  linkDistance = JS("function(d) { return 80 + d.value * 5; }")  
)

# 添加 JavaScript 让图形居中
network_plot <- htmlwidgets::onRender(network_plot, "
  function(el, x) {
    d3.forceSimulation()
      .force('center', d3.forceCenter(el.getBoundingClientRect().width / 2, 
                                      el.getBoundingClientRect().height / 2));
    d3.selectAll('.node text')
      .style('fill', '#000000')
  }
")

  return(network_plot)
}
#######
##setwd("D://博士期间//课题组//博士课题//方法//ATAC//我用的代码//FigR")##本地画图
load("downHRG_figR.d.Rdata")
#load("allHRG_figR.d.Rdata")
#load("upHRG_logFC0.5_figR.d.RData")
#load("ISG_figR.d.Rdata")
#load("downHRG_loFC0.5_figR.d.RData")
label <- read.table("node_label_downHRG.txt",header=F,sep="\t")###去掉密集聚集点中的部分基因
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
####ISG
plotfigRNetwork.2.part(figR.d,
                                              score.cut = 1.5,
                                              weight.edges = TRUE,
                                              legend=FALSE,##TRUE显示图例
                                              charge = -13,
                                              fontSize= 10,
                                             size_dorc = 4,
                                             size_TF = 16,
                                             genes.to.label = node_label_downHRG)
###plotfigRNetwork.2.part.manydorc与plotfigRNetwork.2.part的区别是dorc太多，node太多，原来的node白色描边太粗，所以变细了node的白色边
plotfigRNetwork.2.part.manydorc(figR.d,
                       score.cut = 2,
                       weight.edges = TRUE,
                       legend=FALSE,##TRUE显示图例
                       charge = -1,
                       fontSize= 12,
                       size_dorc = 1,
                       size_TF = 14,
                       genes.to.label = genes.to.label)
###在Rstudio上面不能直接导出pdf####
#先保存html文件
library(htmlwidgets)
saveWidget(plotfigRNetwork.2.part(figR.d, 
                                  score.cut = 2,
                                  weight.edges = TRUE,
                                  legend=TRUE, ## TRUE 显示图例
                                  charge = -15,
                                  fontSize= 12,
                                  size_dorc = 10,
                                  size_TF = 18,
                                  genes.to.label = node_label_downHRG), 
           "downHRG_logFC0.5_network_plot.html", 
           selfcontained = TRUE) 
####保存
saveWidget(plotfigRNetwork.2.part.manydorc(figR.d,
                                           score.cut = 2,
                                           weight.edges = TRUE,
                                           legend=FALSE,##TRUE显示图例
                                           charge = -1,
                                           fontSize= 12,
                                           size_dorc = 1,
                                           size_TF = 14,
                                           genes.to.label = genes.to.label), 
           "upHRG_logFC0.5_network_plot.html", 
           selfcontained = TRUE) 
####ISG
saveWidget(plotfigRNetwork.2.part.manydorc(figR.d,
                                  score.cut = 1.5,
                                  weight.edges = TRUE,
                                  legend=FALSE,##TRUE显示图例
                                  charge = -13,
                                  fontSize= 10,
                                  size_dorc = 4,
                                  size_TF = 16,
                                  genes.to.label = node_label_downHRG), 
           "ISG_network_plot.描边变细.html", 
           selfcontained = TRUE)
######只画个别基因的networkPlot######
figR.d <- figR.d[figR.d$Motif %in% c("KLF10","LEF1","TCF7"),]
library(htmlwidgets)
saveWidget(plotfigRNetwork.2.part(figR.d, 
                                  score.cut = 2,
                                  weight.edges = TRUE,
                                  legend=TRUE, ## TRUE 显示图例
                                  charge = -500,
                                  fontSize= 12,
                                  size_dorc = 14,
                                  size_TF = 18,
                                  genes.to.label = node_label_downHRG), 
           "3TF-downHRG_logFC0.5_network_plot.html", 
           selfcontained = TRUE)
