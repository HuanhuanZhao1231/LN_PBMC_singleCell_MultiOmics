##############
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
# Color Maps
scriptPath_color <- "~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts"
allcolour_atac <- readRDS(paste0(scriptPath_color, "/allcolour_atac.rds")) %>% unlist()
allcolour_RNA <- readRDS(paste0(scriptPath_color, "/allcolour_RNA.rds")) %>% unlist()
broadClustCmap <- readRDS(paste0(scriptPath_color, "/broadClustCmap.rds")) %>% unlist()
##########
atac_proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-proj_N9_MonoSub_ProjHeme5")

############overlapcutoff####
proj_overlap05 <- addPermPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI",
    corCutOff = 0.5,
    overlapCutoff = 0.5,
    addPermutedPval = FALSE,
    seed = 1,
    threads = 16
)

proj_overlap03 <- addPermPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI",
    corCutOff = 0.5,
    overlapCutoff = 0.3,
    addPermutedPval = FALSE,
    seed = 1,
    threads = 16
)

p2g08 <- metadata(atac_proj@peakSet)$Peak2GeneLinks


p2g05 <- metadata(proj_overlap05@peakSet)$Peak2GeneLinks


p2g03 <- metadata(proj_overlap03@peakSet)$Peak2GeneLinks


get_sig <- function(x){

    x <- x[
        !is.na(x$Correlation) &
        x$Correlation>0.5,
    ]

    return(x)
}


sig08 <- get_sig(p2g08)
sig05 <- get_sig(p2g05)
sig03 <- get_sig(p2g03)

nrow(sig08)
nrow(sig05)
nrow(sig03)
[1] 55743
[1] 57423
[1] 57113

gene08 <- unique(sig08$idxRNA)
gene05 <- unique(sig05$idxRNA)
gene03 <- unique(sig03$idxRNA)

length(intersect(gene08,gene05))/length(gene08)
[1] 0.9419114

length(intersect(gene08,gene03))/length(gene08)
[1] 0.939153

#################
perm_oc8_proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-Perm_proj_N9_MonoSub_ProjHeme5_before_P2G")
perm_oc5_proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-Perm_OC5_proj_N9_MonoSub_ProjHeme5_before_P2G")
perm_oc3_proj <- loadArchRProject("~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Save-Perm_OC3_proj_N9_MonoSub_ProjHeme5_before_P2G")

############################################################
## Initialize
############################################################
subclustered_projects <- c(
"T/NK",
"Myeloid",
"Bcells"
)
Perm_sub_dirs <- list(

T/NK =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_T/NK",

Myeloid =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_Myeloid",

Bcells =
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/realCellsubset/Perm_peak2Glinks_Bcells"

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
allPeaksGR <- getPeakSet(perm_oc8_proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["pbmc"]] <- getPeak2GeneLinks_perm(perm_oc8_proj, corCutOff=corrCutoff, FDRCutOff=Inf,PermFDRCutOff=0.05,resolution = 100)[[1]]
coaccessibility_list <- list()
coAccPeaks <- getCoAccessibility(perm_oc8_proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
coAccPeaks$source <- "pbmc"
coaccessibility_list[["pbmc"]] <- coAccPeaks
peak2gene_list <- list()
p2gGR <- getP2G_GR_perm(perm_oc8_proj, corrCutoff=NULL, PermFDRCutOff=NULL,varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
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
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_p2gGR_permOC3.rds"
)


saveRDS(
plot_loop_list,
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_plot_loops_permOC3.rds"
)


saveRDS(
full_coaccessibility,
"~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/scATAC_preprocessing/Perm_results/p2G_allpeak_multilevel_coaccessibility_permOC3.rds"
)
#########
##########################################################################################
# Filter redundant peak to gene links
##########################################################################################

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(perm_oc8_proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# Collapse redundant p2gLinks:
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)
metadata(new_p2g_DF) <- p2gMeta
metadata(perm_oc8_proj@peakSet)$Peak2GeneLinks <- new_p2g_DF
#######上面对3个对象分别运行3次#########
p2gGR_OC8 <- getP2G_GR_perm(perm_oc8_proj, corrCutoff=corrCutoff,PermFDRCutOff=PermFDRCutOff)
#87,686
p2gGR_OC5 <- getP2G_GR_perm(perm_oc5_proj, corrCutoff=corrCutoff,PermFDRCutOff=PermFDRCutOff)
#76,712
p2gGR_OC3 <- getP2G_GR_perm(perm_oc3_proj, corrCutoff=corrCutoff,PermFDRCutOff=PermFDRCutOff)
#77,735
gene08 <- unique(p2gGR_OC8$idxRNA)
gene05 <- unique(p2gGR_OC5$idxRNA)
gene03 <- unique(p2gGR_OC3$idxRNA)

length(intersect(gene08,gene05))/length(gene08)
[1] 0.898039

length(intersect(gene08,gene03))/length(gene08)
[1] 0.907805
> length(intersect(gene08,gene05))
[1] 8825
length(intersect(gene08,gene03))
[1] 8921
length(gene08)
[1] 9827
> length(gene05)
[1] 9091
> length(gene03)
[1] 9258
##################对结果进行绘图展示##############
##############加载包和数据#########
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)


############################################################
## Input
############################################################

p2g08 <- p2gGR_OC8
p2g05 <- p2gGR_OC5
p2g03 <- p2gGR_OC3


############################################################
## 保留有效link
############################################################

p2g08_df <- as.data.frame(p2g08)
p2g05_df <- as.data.frame(p2g05)
p2g03_df <- as.data.frame(p2g03)


p2g08_df <- p2g08_df %>%
    filter(!is.na(Correlation))

p2g05_df <- p2g05_df %>%
    filter(!is.na(Correlation))

p2g03_df <- p2g03_df %>%
    filter(!is.na(Correlation))


head(p2g08_df)

############################################################
############################################################
## merge P2G links
############################################################

merge08_05 <- merge(
    p2g08_df[,c(
        "idxATAC",
        "idxRNA",
        "Correlation"
    )],
    p2g05_df[,c(
        "idxATAC",
        "idxRNA",
        "Correlation"
    )],
    by=c("idxATAC","idxRNA"),
    suffixes=c("_08","_05")
)



merge08_03 <- merge(
    p2g08_df[,c(
        "idxATAC",
        "idxRNA",
        "Correlation"
    )],
    p2g03_df[,c(
        "idxATAC",
        "idxRNA",
        "Correlation"
    )],
    by=c("idxATAC","idxRNA"),
    suffixes=c("_08","_03")
)


nrow(merge08_05)
#72166
nrow(merge08_03)
#71648
#############绘图函数##########
plot_corr <- function(df,x,y,title){


r <- cor(
    df[[x]],
    df[[y]],
    method="pearson"
)


ggplot(
    df,
    aes_string(
        x=x,
        y=y
    )
)+
geom_point(
    size=0.8,
    alpha=0.25
)+
geom_smooth(
    method="lm",
    se=FALSE,
    linewidth=0.8
)+
coord_fixed()+

annotate(
"text",
x=-Inf,
y=Inf,
label=paste0(
"Pearson r = ",
round(r,3)
),
hjust=-0.1,
vjust=1.5,
size=4
)+

labs(
x="Correlation (overlapCutoff=0.8)",
y=y,
title=title
)+

theme_classic()+
theme(
plot.title=element_text(
hjust=0.5
)
)

}
#################
p_corr05 <- plot_corr(
    merge08_05,
    "Correlation_08",
    "Correlation_05",
    "0.8 vs 0.5"
)


p_corr03 <- plot_corr(
    merge08_03,
    "Correlation_08",
    "Correlation_03",
    "0.8 vs 0.3"
)


p_corr <- p_corr05+p_corr03

ggsave(file.path(plotDir,"overlap358_sensitivity_analysis.pdf"),p_corr,width=10,height=5)
#############Figure2-P2G link number
link_num <- data.frame(
    overlapCutoff=c(
        "0.8",
        "0.5",
        "0.3"
    ),
    Links=c(
        length(p2g08_df$idxRNA),
        length(p2g05_df$idxRNA),
        length(p2g03_df$idxRNA)
    )
)


link_num

p_links <- ggplot(
    link_num,
    aes(
        x=overlapCutoff,
        y=Links
    )
)+

geom_col(
    width=0.65
)+

geom_text(
    aes(
        label=scales::comma(Links)
    ),
    vjust=-0.5,
    size=4
)+

labs(
    x="KNN overlap cutoff",
    y="Number of significant peak-gene links"
)+

theme_classic()+
theme(
    text=element_text(size=12)
)


p_links
ggsave(file.path(plotDir,"overlap358_link_number.pdf"),p_links,width=5,height=5)

###########Gene-level reproducibility#####
gene_recovery <- data.frame(

cutoff=c(
    "0.5",
    "0.3"
),

Recovered=c(

length(intersect(
    gene08,
    gene05
))/length(gene08),

length(intersect(
    gene08,
    gene03
))/length(gene08)

)

)


gene_recovery

p_gene <- ggplot(
    gene_recovery,
    aes(
        x=cutoff,
        y=Recovered
    )
)+

geom_col(
    width=0.6
)+

geom_text(
    aes(
        label=paste0(
            round(Recovered*100,1),
            "%"
        )
    ),
    vjust=-0.5,
    size=5
)+

scale_y_continuous(
    limits=c(0,1)
)+

labs(
    x="KNN overlap cutoff",
    y="Gene recovery compared with cutoff=0.8"
)+

theme_classic()


p_gene
ggsave(file.path(plotDir,"overlap358_gene_reproductivity.pdf"),p_gene,width=5,height=5)