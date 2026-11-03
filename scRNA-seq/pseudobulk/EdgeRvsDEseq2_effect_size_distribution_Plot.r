##########effect_size_distribution_edgeR_vs_apeglm

library(tidyverse)

setwd("~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/celltype_DEG/DEG_Plot")
##################################################
# read results
##################################################

edgeR <- read.csv(
"~/scRNA/B/10HC11LN/DEG/data/LNvsHC_EdgeR_DEG_immune.combined.csv"
)


apeglm <- read.csv(
"~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/celltype_DEG/DESeq2_pseudobulk_results/immune.combined_DESeq2_apeglm_DEG.csv"
)

head(edgeR)
head(apeglm)
colnames(edgeR) <- c("gene","logFC","logCPM","PValue","FDR")
lfc_compare <- inner_join(
  
  edgeR %>%
    select(
      gene,
      edgeR_logFC=logFC
    ),
  
  apeglm %>%
    select(
      gene,
      apeglm_logFC=log2FoldChange
    ),
  
  by="gene"
)
lfc_compare <- inner_join(
  
  edgeR %>%
    select(
      gene,
      edgeR_logFC=logFC
    ),
  
  apeglm %>%
    select(
      gene,
      apeglm_logFC=log2FoldChange
    ),
  
  by="gene"
)

lfc_long <- lfc_compare %>%
  
  pivot_longer(
    cols=c(
      edgeR_logFC,
      apeglm_logFC
    ),
    
    names_to="Method",
    values_to="logFC"
  )


head(lfc_long)

###densityPlot#####
library(ggplot2)


p_density <- ggplot(
  
  lfc_long,
  
  aes(
    x=logFC,
    color=Method
  )

)+


geom_density(
  linewidth=1
)+


geom_vline(
  xintercept=c(-0.5,0.5),
  linetype="dashed"
)+


theme_classic(
  base_size=14
)+


labs(
  
  x="Estimated log2 fold change",
  
  y="Density",
  
  color=NULL
  
)



p_density

ggsave(
"Figure_0.5effect_size_distribution_edgeR_vs_apeglm.pdf",
p_density,
width=6,
height=4
)