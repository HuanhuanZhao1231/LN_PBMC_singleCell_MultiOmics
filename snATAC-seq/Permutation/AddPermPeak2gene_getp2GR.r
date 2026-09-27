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
scriptPath <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/perm_functions.R"))
# Set Threads to be used
addArchRThreads(threads = 16)
# set working directory (The directory of the full preprocessed archr project)
plotDir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2gLink_plots"
tableDir <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/supplemental_tables"
# Color Maps
scriptPath_color <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
allcolour_atac <- readRDS(paste0(scriptPath_color, "/allcolour_atac.rds")) %>% unlist()
allcolour_RNA <- readRDS(paste0(scriptPath_color, "/allcolour_RNA.rds")) %>% unlist()
broadClustCmap <- readRDS(paste0(scriptPath_color, "/broadClustCmap.rds")) %>% unlist()
##########
atac_Proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_MonoSub_ProjHeme5")
perm_proj <- addPermPeak2GeneLinks(
    ArchRProj=atac_proj,
    reducedDims="IterativeLSI",
    corCutOff=0.5,
    addPermutedPval=TRUE,
    nperm=1000,
    seed=123,
    threads=16
)

p2g_perm <- metadata(perm_proj@peakSet)$Peak2GeneLinks


summary(p2g_perm$PermFDR)


sum(
p2g_perm$Correlation>0.5 &
p2g_perm$PermFDR<0.05,
na.rm=TRUE
)
setwd("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset")
saveArchRProject(
    ArchRProj = perm_proj,
    outputDirectory="Save-Perm_proj_N9_MonoSub_ProjHeme5_before_P2G",
    load=FALSE
)
#########################################重新加载perm_proj########
perm_proj <- 
loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-Perm_OC5_proj_N9_MonoSub_ProjHeme5_before_P2G")

colnames(
metadata(perm_proj@peakSet)$Peak2GeneLinks
)
subclustered_projects <- c(
"T/NK",
"Myeloid",
"Bcells"
)

################3个subGroup重新addpermPeak2Genelinks##########
run_subproject_permP2G <- function(
    project_dir,
    name,
    save_dir,
    nperm=1000
){


    message(
        "Loading ",
        name
    )


    proj <- loadArchRProject(
        project_dir,
        force=TRUE
    )


    message(
        "Running permutation P2G for ",
        name
    )


    proj <- addPermPeak2GeneLinks(
        ArchRProj = proj,
        reducedDims = "IterativeLSI",
        corCutOff = 0.5,
        addPermutedPval = TRUE,
        nperm = nperm,
        seed = 123,
        threads = 16
    )


    saveArchRProject(
        ArchRProj = proj,
        outputDirectory = save_dir,
        load = FALSE
    )


    return(
        list(
            proj=proj
        )
    )

}
lymph_result <- run_subproject_permP2G(
    
    project_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_T/NK",
    
    name="T/NK",
    
    save_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_T/NK",
    
    nperm=1000
)
myeloid_result <- run_subproject_permP2G(
    
    project_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_Myeloid",
    
    name="Myeloid",
    
    save_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_Myeloid",
    
    nperm=1000
)
Bcells_result <- run_subproject_permP2G(
    
    project_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/peak2Glinks_Bcells",
    
    name="Bcells",
    
    save_dir =
    "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_Bcells",
    
    nperm=1000
)

######################
Perm_sub_dirs <- list(

T/NK =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC5_peak2Glinks_T/NK",

Myeloid =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC5_peak2Glinks_Myeloid",

Bcells =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_OC5_peak2Glinks_Bcells"

)
################重新生成P2G和loop
############################################################
## Initialize
############################################################
subclustered_projects <- c(
"T/NK",
"Myeloid",
"Bcells"
)

peak2gene_list <- list()

plot_loop_list <- list()

coaccessibility_list <- list()



############################################################
## PBMC
############################################################
# P2G definition cutoffs
corrCutoff <- 0.5       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25
PermFDRCutOff <- 0.05
# Coaccessibility cutoffs
coAccCorrCutoff <- 0.4  # Default in getCoAccessibility is 0.5

# Get all peaks
allPeaksGR <- getPeakSet(perm_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["pbmc"]] <- getPeak2GeneLinks_perm(perm_proj, corCutOff=corrCutoff, FDRCutOff=Inf,PermFDRCutOff=0.05,resolution = 100)[[1]]
coaccessibility_list <- list()
coAccPeaks <- getCoAccessibility(perm_proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
coAccPeaks$source <- "pbmc"
coaccessibility_list[["pbmc"]] <- coAccPeaks
peak2gene_list <- list()
p2gGR <- getP2G_GR_perm(perm_proj, corrCutoff=NULL, PermFDRCutOff=NULL,varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
p2gGR$source <- "pbmc"
peak2gene_list[["pbmc"]] <- p2gGR

#########3个subProject###########
for(subgroup in subclustered_projects){


    message(
        "Processing ",
        subgroup
    )


    sub_proj <-
        loadArchRProject(
            Perm_sub_dirs[[subgroup]],
            force=TRUE
        )



    ################################################
    # P2G
    ################################################


    subP2G <-

    getP2G_GR_perm(
        sub_proj,
        corrCutoff=NULL,
        varCutoffATAC=-Inf,
        varCutoffRNA=-Inf,
        PermFDRCutOff=NULL,
        filtNA=FALSE
    )


    subP2G$source <-
        subgroup


    peak2gene_list[[subgroup]] <-
        subP2G



    ################################################
    # loops
    ################################################


    plot_loop_list[[subgroup]] <-

    getPeak2GeneLinks_perm(
        sub_proj,
        corCutOff=0.5,
        FDRCutOff=Inf,
        PermFDRCutOff=0.05,
        resolution=100
    )[[1]]



    ################################################
    # coAccessibility
    ################################################


    coAccPeaks <-

    getCoAccessibility(
        sub_proj,
        corCutOff=0.5,
        returnLoops=TRUE
    )[[1]]


    coAccPeaks$source <-
        subgroup


    coaccessibility_list[[subgroup]] <-
        coAccPeaks

}

###############合并###############
############################################################
## Merge
############################################################


full_p2gGR <-

as(
    peak2gene_list,
    "GRangesList"
) %>%
unlist()



full_coaccessibility <-

as(
    coaccessibility_list,
    "GRangesList"
) %>%
unlist()



############################################################
## fix idxATAC
############################################################


idxATAC <-
peak2gene_list[["pbmc"]]$idxATAC


names(idxATAC) <-
peak2gene_list[["pbmc"]]$peakName



full_p2gGR$idxATAC <-

idxATAC[
    full_p2gGR$peakName
]



saveRDS(
full_p2gGR,
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_p2gGR_permOC8.rds"
)


saveRDS(
plot_loop_list,
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_plot_loops_permOC8.rds"
)


saveRDS(
full_coaccessibility,
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_coaccessibility_permOC8.rds"
)
##这里是将整体peak2gene和每个subgroup的peak2gene合起来的总体p2gGR
full_p2gGR <- readRDS("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_p2gGR_permOC5.rds")
plot_loop_list <- readRDS("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_plot_loops_permOC5.rds")
full_coaccessibility <- readRDS("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_coaccessibility_permOC5.rds")
##########################################################################################
# Upset plot of number of peak to gene links identified per group


library(UpSetR)

groups <- unique(full_p2gGR$source)
sub_p2gGR <- full_p2gGR[!is.na(full_p2gGR$Correlation)]

upset_list <- lapply(groups, function(g){
  gr <- sub_p2gGR[sub_p2gGR$source == g & 
    sub_p2gGR$Correlation > corrCutoff &
    sub_p2gGR$PermFDR < PermFDRCutOff & 
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

pdf(paste0(plotDir, "/P2g_upset_plot_linked_peaks_permOC5.pdf"), width=10, height=5)
plotUpset(upset_list, main.bar.color="royalblue1", keep.order=FALSE)
dev.off()

##########################################################################################
# Filter redundant peak to gene links
##########################################################################################

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(perm_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)
metadata(new_p2g_DF) <- p2gMeta
metadata(perm_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF
# Plot some comparisons between linked peaks and unlinked peaks
##########################################################################################
p2gGR <- getP2G_GR_perm(perm_proj, corrCutoff=corrCutoff,PermFDRCutOff=PermFDRCutOff)
#PermFDR0.05_corrCutoff_87686P2G
#No_Perm结果_91479P2G,corrCutoff0.5
#No_Perm结果_总peaks数目#242637
all_linked_peaks <- p2gGR$peakName %>% unique()
#No_perm_53853
#Perm_52260
unlinked_peaks <- allPeaksGR[allPeaksGR$peakName %ni% all_linked_peaks]$peakName
#188784
#Perm_190377
#Most peaks (188784, 77.80512%) were not linked to any gene, consistent with the expected small effect size of most CREs32.
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
pdf(paste0(plotDir, "/linked_vs_unlinked_GC_content_p2G.pdf"), width=6, height=6)
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
pdf(paste0(plotDir, "/Perm5_linked_vs_unlinked_distance_to_TSS_p2G.pdf"), width=6, height=5)
p
dev.off()
##########################################################################################
# Conservation of linked peaks vs unlinked peaks in each subgrouping
##########################################################################################

library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19

filt_full_p2gGR <- full_p2gGR[full_p2gGR$Correlation > corrCutoff & 
    full_p2gGR$VarQATAC > varCutoffATAC & 
    full_p2gGR$VarQRNA > varCutoffRNA]

filt_full_p2gGR$p2gName <- paste0(filt_full_p2gGR$peakName, "_", filt_full_p2gGR$symbol)

# For each source, get the phastCons100 way median conservation in linked peaks
allPeaksCons <- gscores(phast, allPeaksGR, summaryFun=mean) # Mean is the default summaryFun

# Get conservation of linked and unlinked peaks from each dataset
p2g_groups <- unique(filt_p2gGR$source)

cons_df <- lapply(p2g_groups, function(g){
  sub_p2g_names <- filt_full_p2gGR$peakName[filt_full_p2gGR$source == g] %>% unique()
  lcons <- allPeaksCons$default[allPeaksCons$peakName %in% sub_p2g_names]
  data.frame(
    group=rep(g, times=length(lcons)), 
    conservation=lcons
    )
  }) %>% do.call(rbind,.)

# Add conservation of peaks with no link
ulcons <- allPeaksCons$default[allPeaksCons$peakName %ni% filt_full_p2gGR$peakName]
cons_df <- rbind(cons_df, data.frame(
  group=rep("unlinked", times=length(ulcons)),
  conservation=ulcons
  ))

cons_pvals <- lapply(p2g_groups, function(g){
  glcons <- cons_df[cons_df$group == g, "conservation"]
  ulcons <- cons_df[cons_df$group == "unlinked", "conservation"]
  c(linked_mean=mean(glcons, na.rm=TRUE), unlinked_mean=mean(ulcons, na.rm=TRUE), pval=wilcox.test(glcons, ulcons)$p.value)
  }) %>% do.call(rbind,.) %>% as.data.frame()
rownames(cons_pvals) <- p2g_groups

cons_cmap <- rep("royalblue1", times=length(p2g_groups))
names(cons_cmap) <- p2g_groups
cons_cmap <- c(cons_cmap, unlinked="grey")

# Order by decreasing mean conservation
cons_order <- rownames(cons_pvals[order(cons_pvals$linked_mean, decreasing=TRUE),])
cons_df$group <- factor(cons_df$group, levels=c(cons_order, "unlinked"), ordered=TRUE)

p <- (
  ggplot(cons_df, aes(x=group, y=conservation, fill=group), color="black")
  + geom_boxplot(alpha=1.0)
  + scale_y_continuous(limits=c(0.0,1.05), expand=c(0,0))
  + scale_color_manual(values=cons_cmap)
  + scale_fill_manual(values=cons_cmap)
  + xlab("")
  + ylab("phastCons.100")
  + theme_BOR(border=FALSE)
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          legend.position="none", # Remove legend
          axis.text.x = element_text(angle=90, hjust=1)) 
)

pdf(paste0(plotDir, "/PermOC5_linked_vs_unlinked_peak_conservationp2G.pdf"), width=6, height=5)
p
dev.off()

##########################################################################################

# Identify 'highly regulated' genes and 'highly-regulating' peaks
##########################################################################################
############################################################
## 2. Count linked peaks per gene
############################################################
p2gGR <- getP2G_GR_perm(perm_proj, corrCutoff=corrCutoff,PermFDRCutOff=PermFDRCutOff)
p2g_sig <- p2gGR
library(dplyr)
gene_peak_count <- as.data.frame(p2g_sig) %>%
    group_by(idxRNA) %>%
    summarise(
        nLinkedPeak = n_distinct(idxATAC)
    ) %>%
    arrange(desc(nLinkedPeak))
head(gene_peak_count)
################绘制

addArchRThreads(threads = 1)

load("~/scRNA/B/10HC11LN/3score_labled15.8immune.combinedCelltypeGroup.RData")  
# Get all expressed genes:
count.mat <- Seurat::GetAssayData(object=immune.combined, slot="counts")
minUMIs <- 1
minCells <- 2
valid.genes <- rownames(count.mat[rowSums(count.mat > minUMIs) > minCells,])

# Get distribution of peaks to gene linkages and identify 'highly-regulated' genes
p2gGR <- getP2G_GR_perm(perm_proj, corrCutoff=corrCutoff, PermFDRCutOff=PermFDRCutOff)
p2gFreqs <- getFreqs(p2gGR$symbol)
valid.genes <- c(valid.genes, unique(p2gGR$symbol)) %>% unique()

noLinks <- valid.genes[valid.genes %ni% names(p2gFreqs)]
zilch <- rep(0, length(noLinks))
names(zilch) <- noLinks
p2gFreqs <- c(p2gFreqs, zilch)
x <- 1:length(p2gFreqs)
rank_df <- data.frame(npeaks=p2gFreqs, rank=x)
write.table(rank_df,file.path(tableDir,"Perm_OC5_Peak2Gene_rank_df.txt"),quote=FALSE,row.names=TRUE,col.names=TRUE) 
## 4. Identify elbow point
############################################################
#因为你的gene peak数量通常长尾分布：后面的长尾会影响elbow。建议：只考虑：npeaks >= 5
rank_df_5 <- rank_df %>%
    filter(npeaks>=5)
###计算elbow#
## Elbow detection by maximum distance method
#########################################################
elbow_df <- rank_df_5
# normalize coordinates
x <- elbow_df$rank
y <- elbow_df$npeaks
x_norm <- (x-min(x))/(max(x)-min(x))
y_norm <- (y-min(y))/(max(y)-min(y))
# line connecting first and last points
x1 <- x_norm[1]
y1 <- y_norm[1]

x2 <- x_norm[length(x_norm)]
y2 <- y_norm[length(y_norm)]
# distance from each point to line
distance <- abs(
    (y2-y1)*x_norm -
    (x2-x1)*y_norm +
    x2*y1 -
    y2*x1
) /
sqrt(
    (y2-y1)^2+
    (x2-x1)^2
)
elbow_index <- which.max(distance)
elbow_index
#Step 2. 查看对应peak cutoff
elbow_cutoff <- rank_df$npeaks[elbow_index]
elbow_cutoff
#25
#########画elbowPlot先看一下###
library(ggplot2)
p <- ggplot(
    rank_df_5,
    aes(
        x=rank,
        y=npeaks
    )
)+
geom_line()+
geom_point(
    size=0.8
)+
geom_point(
    data=rank_df_5[elbow_index,],
    aes(
        x=rank,
        y=npeaks
    ),
    size=4
)+
geom_hline(
    yintercept=elbow_cutoff,
    linetype="dashed"
)+
theme_classic()+
labs(
    x="Gene rank",
    y="Number of significant peak-to-gene links",
    title=paste0(
        "Elbow cutoff = ",
        elbow_cutoff,
        " peaks"
    )
)
ggsave("ran_df_5_elbowPlot.pdf",p)
#Step 5. Threshold sensitivity（Reviewer明确要求）
#计算：
10,15,20,30。
############################################################
## 5. HRG sensitivity analysis
############################################################
thresholds <- c(10,15,20,25,30)
HRG_summary <- data.frame(
    cutoff=thresholds,
    nHRG=sapply(
        thresholds,
        function(x){
            sum(
                rank_df$npeaks >= x
            )
        }
    )
)

HRG_summary
#> HRG_summary
#  cutoff nHRG
#1     10 2754
#2     15 1682
#3     20 1040
#4     25  622
#5     30  386
#######6.HRG稳定性
###################################################
## Jaccard similarity between HRG sets
###################################################

cutoffs <- c(10,15,20,25,30)


jaccard_mat <- matrix(
    NA,
    nrow=length(cutoffs),
    ncol=length(cutoffs)
)

rownames(jaccard_mat)<-
colnames(jaccard_mat)<-
paste0("HRG",cutoffs)


for(i in seq_along(cutoffs)){
    
    for(j in seq_along(cutoffs)){
        
        A <- HRG_list[[i]]
        B <- HRG_list[[j]]
        
        jaccard_mat[i,j] <-
            length(intersect(A,B))/
            length(union(A,B))
    }
}


round(jaccard_mat,3)
> round(jaccard_mat,3)
      HRG10 HRG15 HRG20 HRG25 HRG30
HRG10 1.000 0.611 0.378 0.226 0.140
HRG15 0.611 1.000 0.618 0.370 0.229
HRG20 0.378 0.618 1.000 0.598 0.371
HRG25 0.226 0.370 0.598 1.000 0.621
HRG30 0.140 0.229 0.371 0.621 1.000
# Cutoff for defining highly regulated genes (HRGs) determined using elbow rule
cutoff <- 25
# Save HRGs as table
hrg_df <- rank_df_5[rank_df_5$npeaks >= cutoff,]
hrg_df$gene <- rownames(hrg_df)
write.table(hrg_df, file=paste0(tableDir, "/Perm_OC5_cutoff25_p2G_HRG_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 
 
#################################正式画带gene label的elbowPlot########
# Rank plot:(拒绝长尾效应，选用npeaks>5的gene作图)

rank_df_5$color <- ifelse(rank_df_5$npeaks >= cutoff, "royalblue1", "black")

# Label select super enhancer - associated genes
label_genes <- c(
"JDP2","ZEB2","CEBPB","FOS","BCL11B","BACH2","TCF7",
"RUNX3","ETS1", "IFI30","ETS2","IL6R","IFNGR2",
"KLF4","IRF8","BLK","PAX5","SPI1","PRDM1"
)

rank_df_5$label <- ifelse(rownames(rank_df_5) %in% label_genes, rownames(rank_df_5), "")
p <- (ggplot(data=rank_df_5, aes(x=rank, y=npeaks))
+ geom_line()+
geom_point(
    size=0.5
)+
geom_point(
    data=rank_df_5[elbow_index,],
    aes(
        x=rank,
        y=npeaks
    ),
    size=4
)+
geom_hline(
    yintercept=elbow_cutoff,
    linetype="dashed"
)
  + geom_point_rast(aes(color=color))
  + ggrepel::geom_text_repel(
          data = rank_df_5[rank_df_5$label != "",], aes(x=rank, y=npeaks, label=label), 
          size = 3,
          nudge_x = 2,
          #direction = "x",
          hjust = "outward",
          segment.size = 0.1,
          box.padding=0.5,
          min.segment.length = 0, # draw all segments
          max.overlaps = Inf, # draw all labels
          color = "black")
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

pdf(paste0(plotDir, "/Perm_OC5_cutoff25_nLinkedPeaksPerGene_rastr_p2G点在上.pdf"))
p
dev.off()
# Plot barplot of how many linked peaks per gene
thresh <- 30
threshNpeaks <- rank_df$npeaks
threshNpeaks[threshNpeaks>thresh] <- thresh
nLinkedPeaks <- getFreqs(threshNpeaks)

df <- data.frame(nLinkedPeaks=as.integer(names(nLinkedPeaks)), nGenes=nLinkedPeaks)

pdf(paste0(plotDir, "/p2G_nGenes_with_nLinkedPeaks.pdf"), width=8, height=6)
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

pdf(paste0(plotDir, "/P2G_nPeaks_with_nLinkedGenes.pdf"), width=5, height=6)
qcBarPlot(df, cmap="royalblue1", barwidth=0.9, border_color=NA) + geom_vline(xintercept=median(g2p_rank_df$npeaks), linetype="dashed")
dev.off()

##########################################################################################

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

pdf(paste0(plotDir, "/peak_types_barplot_p2G.pdf"), width=5, height=6)
stackedBarPlot(melt_df, xcol=2, fillcol=1, ycol=3, 
  cmap=getColorMap(cmaps_BOR$sambaNight, n=4, type="quantitative"), barwidth=0.9)
dev.off()

##########################################################################################
#看peak2gene linakage的linkedpeak链接的gene,是nearest Gene的比例##
df <- data.frame(
  peakName=all_linked_peaks, 
  label=c(rep("linked", length(all_linked_peaks)))
)
df$nearestGene <- allPeaksGR[df$peakName]$nearestGene
names(p2gGR) <- mcols(p2gGR)$peakName
df$symbol <- p2gGR[df$peakName]$symbol
################添加新列
df <- df %>%
  mutate(match = ifelse(nearestGene == symbol, 1, 0))
#将df$match 列中的 NA 值替换为 0####后来发现allPeaksGR[df$peakName]$nearestGene有50多个是NA,暂认为不是nearestGene,所以赋值0
df$match[is.na(df$match)] <- 0
# 计算每个值的数量
match_counts <- table(df$match)
# 计算每个值的比例
match_proportions <- prop.table(match_counts)
# 查看结果
match_proportions
# 将比例转换为数据框
match_proportions_df <- as.data.frame(match_proportions)
colnames(match_proportions_df) <- c("match", "proportion")
pdf("linked_peaks_nearestGene_symbol_propotion.pdf",height=4,width=4)
# 使用 ggplot2 创建条形图
ggplot(match_proportions_df, aes(x = as.factor(match), y = proportion, fill = as.factor(match))) +
  geom_bar(stat = "identity",width=0.8) +
  scale_fill_manual(values = c("0" = "#69B1D6", "1" = "#D1938A"), labels = c("0" = "Not Match", "1" = "Match")) +
  labs(title = "Proportion of Match vs Not Match",
       x = "Match Status",
       y = "Proportion",
       fill = "Match Status") +
  theme_classic()
###使用 ggplot2 创建饼图
pdf(file.path(plotDir,"linked_peaks_nearestGene_symbol_pieplot.pdf"),height=4,width=4)
match_proportions_df$label <- paste0(match_proportions_df$match, ": ", round(match_proportions_df$proportion * 100, 1), "%")
ggplot(match_proportions_df, aes(x = "", y = proportion, fill = as.factor(match))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("0" = "#69B1D6", "1" = "#D1938A"), labels = c("0" = "Not Match", "1" = "Match")) +
  labs(title = "Proportion of Match vs Not Match",
       fill = "Match Status") +
  theme_void()+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))
dev.off()
##########################################################################################
##########################################################################################
# Plot Peak2Gene heatmap
##########################################################################################
################################
for (nclust in c(7,8,9,11,12)){
#nclust <- 8
p <- plotPeak2GeneHeatmap(
  perm_proj, 
  corCutOff = corrCutoff, 
  groupBy="FineClust", 
  nPlot = 1000000, returnMatrices=FALSE, 
  k=nclust, seed=1, palGroup=allcolour_atac
  )
pdf(paste0(plotDir, sprintf("/peakToGeneHeatmap_LabelFineClust_k%s_p2G.pdf",nclust)), width=16, height=10)
print(p)
dev.off()
################################
# Need to force it to plot all peaks if you want to match the labeling when you 'returnMatrices'.
p2gMat <- plotPeak2GeneHeatmap(
  perm_proj, 
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
matches <- getMatches(perm_proj, "Motif")
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

pdf(paste0(plotDir, sprintf("/enrichedMotifs_k%s_p2gHM_p2G.pdf",nclust)), width=12, height=12)
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
outfile <- paste0(plotDir, sprintf("/p2G_finaltop250_genes_kclust_k%s.tsv", nclust))
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
outfile <- paste0(plotDir, sprintf("/p2G_finaltop250_genes_GOresults_kclust_k%s.tsv", nclust))
names(GOresults) <- paste0("cluster_", kclust)
write.table(GOresults, file=outfile, quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)
# Plots of GO term enrichments:
pdf(paste0(plotDir, sprintf("/p2G_kclust_GO_3termsBPonlyBarLim_k%s_250gene.pdf", nclust)), width=10, height=2.5)
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
saveRDS(young_se_dt,file.path(tableDir,"cutoff25_hrg_SE_enrichment_young_se.rds"))
# Convert refseq IDs to gene symbols
library(org.Hs.eg.db)
library(biomaRt)
#biomart需要联网，所以在登陆节点将youngse的refseq转成symbol保存refseq_to_symbol
refseq_to_symbol <- read.table("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/supplemental_tables/young_se_biomart_refseq_to_symbol.txt",header=1)
#mart <- useMart("ensembl","hsapiens_gene_ensembl")###需要联网，在登陆节点进行
#refseq_to_symbol <- biomaRt::select(
#  org.Hs.eg.db, 
#  keys=unique(young_se_dt$refseq), 
#  columns=c("REFSEQ", "SYMBOL"), 
#  keytype="REFSEQ"
#  )##biomaRt::select要定义biomaRt的select功能，因为dplyr等其他函数也有这个function,不定义可能默认其他函数的select,会报错
#biomart需要联网，所以在登陆节点将youngse的refseq转成symbol保存refseq_to_symbol
#Error in UseMethod("select") :
#  no applicable method for 'select' applied to an object of class "c('OrgDb', 'AnnotationDb', 'envRefClass', '.environment', 'refClass', 'environment', 'refObject', 'AssayData')"
ref_to_sym <- refseq_to_symbol$SYMBOL
names(ref_to_sym) <- refseq_to_symbol$REFSEQ
young_se_dt$symbol <- ref_to_sym[young_se_dt$refseq]
young_se_dt <- young_se_dt[!is.na(young_se_dt$symbol)]

# Adams super enhancers 2015 (Mouse hair follicles)
fuchs_se_files <- list.files(
  path="/oak/stanford/groups/wjg/boberrey/hairATAC/analyses/resources/superEnhancers/Adams2015Nature_SEs",
  pattern="*.csv$",
  full.names = TRUE
)
cell_source <- str_replace(basename(fuchs_se_files), "\\.csv$", "")
names(fuchs_se_files) <- cell_source

fuchs_se_dt <- lapply(names(fuchs_se_files), function(x){
  dt <- fread(fuchs_se_files[x], header=FALSE, sep=",", skip=1)
  dt$source <- x
  dt
  }) %>% rbindlist()
colnames(fuchs_se_dt) <- c("chr", "start", "end", "enhRank", "mouse_gene", "overlaps_HFSCs", "source")

mouse_to_human <- convertMouseGeneList(fuchs_se_dt$mouse_gene)
convert_M2H <- mouse_to_human$HGNC.symbol
names(convert_M2H) <- mouse_to_human$MGI.symbol
fuchs_se_dt$symbol <- convert_M2H[fuchs_se_dt$mouse_gene]
fuchs_se_dt <- fuchs_se_dt[!is.na(symbol)]

# Get the top N super enhancers from each source:
topN <- 250

se_sources <- unique(young_se_dt$source)
top_SEs_young <- lapply(se_sources, function(s){
  young_se_dt[source == s] %>% arrange(enhRank) %>% slice(1:topN) %>% pull(symbol)
  })
names(top_SEs_young) <- se_sources

#se_sources <- unique(fuchs_se_dt$source)
#top_SEs_fuchs <- lapply(se_sources, function(s){
#  fuchs_se_dt[source == s] %>% arrange(enhRank) %>% slice(1:topN) %>% pull(symbol)
#  })
#names(top_SEs_fuchs) <- se_sources

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
##源代码将pval转换为mlog10padj再画图，现在要求用-log10pval画图，所以重新画一下
plot_df_pval <- data.frame(rank=1:length(pvals), pval=pvals)
rownames(plot_df_pval) <- names(pvals)
write.table(plot_df_pval, file=paste0(tableDir, "/Perm_OC5_cutoff_25_p2G_hrg_young_SE_enrichment_pval_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 
# Plot enrichment of SE's from different sources:
plot_df <- data.frame(rank=1:length(pvals), mlog10padj=-log10(p.adjust(10**-pvals, method="fdr")))
rownames(plot_df) <- names(pvals)
# Save table
write.table(plot_df, file=paste0(tableDir, "/Perm_OC5_cutoff_25_p2G_hrg_young_SE_enrichment_mlog10padj_table.tsv"), quote=FALSE, sep="\t", col.names=NA, row.names=TRUE) 

to_label <- c("CD20","BI_CD4p_CD25-_Il17-_PMAstim_Th", "BI_CD4_Naive_Primary_8pool", "CD19_primary","CD14", "CD3",
 "UCSD_Lung", "BI_Brain_Hippocampus_Middle", "HeLa", "BI_Adipose_Nuclei","UCSD_Spleen", "CD34_fetal", "HepG2")
plot_df_pval$label <- ifelse(rownames(plot_df_pval) %in% to_label, rownames(plot_df_pval), "")
plot_df_pval$color <- "royalblue1"

p <- (ggplot(data=plot_df_pval, aes(x=rank, y=pval))###原始代码用data=plot_df, y=mlog10padj
  + ggrepel::geom_text_repel(
          data = plot_df_pval[plot_df_pval$label != "",], aes(x=rank, y=pval, label=label), 
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
  + ylab("Hypergeometric Enrichment -log10(Pval)")
  + xlab("")
  + theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
          aspect.ratio = 1.0,
          legend.position = "none", # Remove legend
          axis.text.x = element_text(angle = 90, hjust = 1))
)

pdf(paste0(plotDir, "/Perm_OC5_pval_cutoff_25_youngSE_enrichment_by_source_P2G.pdf"))
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


###################################################################################################
# Plot Tracks of ALL peak to gene links for select super-enhancer genes
###################################################################################################

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

# Tracks of genes:
# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(perm_proj))

mPromoterGR <- promoterGR[promoterGR$symbol %in% label_genes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% label_genes]

# Restrict to only loops linking genes of interest
plotLoops <- getPeak2GeneLinks_perm(perm_proj, corCutOff=corrCutoff, PermFDRCutoff= PermFDRCutOff,resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Bracket plot regions around SNPs
plotRegions <- lapply(label_genes, function(x){
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  gr <- GRanges(
      seqnames = seqnames(gr)[1],
      ranges = IRanges(start=lims[1], end=lims[2])
    )
  gr
  }) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
  width=width(plotRegions) + 0.05*width(plotRegions), 
  fix="center")


# Tracks of genes:
p <- plotBrowserTrack(
    ArchRProj = perm_proj, 
    groupBy = "LNamedClust", 
    useGroups = unlist(atac.NamedClust)[atacOrder],
    pal = atacLabelClustCmap,
    plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"), # Doesn't change order...
    sizes = c(7, 0.2, 1.25, 2.5),
    geneSymbol = label_genes, 
    region = plotRegions, 
    loops = plotLoops,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Super-Enhancers.pdf", 
    ArchRProj = perm_proj, 
    addDOC = FALSE, 
    width = 6, height = 7)


##########################################################################################
# Violin plots of (integrated) RNA expression for select genes
##########################################################################################

# WARNING: this may not work if you have already assigned the new p2glinks to the project
GImat <- getMatrixFromProject(perm_proj, useMatrix="GeneIntegrationMatrix")
data_mat <- assays(GImat)[[1]]
rownames(data_mat) <- rowData(GImat)$name
sub_mat <- data_mat[label_genes,]

# These DO NOT match the order of the above matrix by default
grouping_data <- data.frame(cluster=factor(perm_proj$LNamedClust, 
  ordered=TRUE, levels=unlist(atac.NamedClust)[atacOrder]))
rownames(grouping_data) <- getCellNames(perm_proj)
sub_mat <- sub_mat[,rownames(grouping_data)]

dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)

pList <- list()
for(gn in label_genes){
  df <- data.frame(grouping_data, gene=sub_mat[gn,])
  # Sample to no more than 500 cells per cluster
  df <- df %>% group_by(cluster) %>% dplyr::slice(sample(min(500, n()))) %>% ungroup()
  df <- df[df$cluster %in% unlist(atac.NamedClust)[atacOrder],]

  covarLabel <- "cluster"  

  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=cluster, y=gene, fill=cluster))
    + geom_violin(aes(fill=cluster), adjust = 1.0, scale='width', position=dodge)
    #+ geom_jitter(aes(group=Sample), size=0.025, 
    #  position=position_jitterdodge(seed=1, jitter.width=0.05, jitter.height=0.0, dodge.width=dodge_width))
    #+ stat_summary(fun="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
    # width=0.75, position=dodge,show.legend = FALSE)
    + scale_color_manual(values=atacLabelClustCmap, limits=names(atacLabelClustCmap), name=covarLabel, na.value="grey")
    + scale_fill_manual(values=atacLabelClustCmap)
    + guides(fill=guide_legend(title=covarLabel), 
      colour=guide_legend(override.aes = list(size=5)))
    + ggtitle(gn)
    + xlab("")
    + ylab("Integrated RNA Expression")
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  pList[[gn]] <- p
}

pdf(paste0(plotDir, "/Expression_Violin_byClust.pdf"), width=10, height=4)
pList
dev.off() 
