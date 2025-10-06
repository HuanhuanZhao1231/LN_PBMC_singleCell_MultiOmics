## Load libraries
```{r}
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(edgeR)
library(limma)
#library(G0.db)
```
## Reading in the data

```{r}
in_dir = "~/scRNA-seq/B/10HC11LN/sample_sum_counts.csv"
cts <- read.csv(in_dir, row.names="X")

```

## preparing a DGEList object

```{r}

## read in table with column data
#coldata <- read.csv(paste0(out_dir, "EdgeR_coldata_for_samplesums.csv"), row.names = 1)
clinical_data <- read.csv("~/scRNA-seq/B/10HC11LN/clinical_data.csv")

obj <- DGEList(counts=cts, group=clinical_data$IE, samples = clinical_data)
obj$samples
head(obj$counts)

#filter out lowly expressed genes
keep <- filterByExpr(obj, group = clinical_data$IE,
             min.count = 30, min.total.count = 300, large.n = 4, min.prop = 0.6)
table(keep)
obj <- obj[keep, , keep.lib.sizes=FALSE]

## Normalisation for RNA composition using TMM (trimmed mean of M-values)
obj <- calcNormFactors(obj)
obj$samples
```

## Setting up the model

```{r}
## Set up the design matrix
Sample <- factor(clinical_data$Patient_ID)
IE <- factor(clinical_data$IE)
IE
IE <- relevel(IE, "IE2")

design <- model.matrix(~IE)

## Estimate dispersion (estimates common dispersion, trended dispersions and tagwise dispersions in one run)
obj <- estimateDisp(obj, design)
plotBCV(obj)
```
## Calculate Differential Expression

```{r}
# exact test (only for single-factor experiments)
et <- exactTest(obj)
topTags(et)

# Number of up/downregulated genes at 5% FDR
summary(decideTests(et))

plotMD(et)
abline(h=c(-1,1), col ="blue")
```


## Exporting results
```{r}
out_dir = "~/scRNA-seq/B/10HC11LN/output/pseudobulk/"
write.csv(as.data.frame(topTags(et, n=Inf)), file=paste0(out_dir, "LNvsHC_EdgeR_samplesums_filtered.csv"))
```

## Convert to Entrez Gene IDs (NCBI IDs)

```{r}
## Prepare Mart and Conversion Table

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_list = unique(row.names(cts))
conversion_table=getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = gene_list, bmHeader = T, mart = mart)
colnames(conversion_table)=c("hgnc_symbol","ncbi_gene_id")
#conversion_table$ncbi_gene_id = ifelse(conversion_table$ncbi_gene_id == "", conversion_table$hgnc_symbol, conversion_table$ncbi_gene_id)

## Exchange columns
setDT(cts, keep.rownames = 'hgnc_symbol')
cts_new = merge.data.frame(cts, conversion_table, by = "hgnc_symbol")
cts_new$hgnc_symbol = NULL
cts_new = cts_new[c(15, 1:14)]

# Find and remove missing values
NAcount = is.na(cts_new$ncbi_gene_id)
table(NAcount)["TRUE"]
cts_complete <- na.omit(cts_new, cols="ncbi_gene_id")

# Find and remove duplicates
duplicateCount = cts_complete[duplicated(cts_complete$ncbi_gene_id),]
View(duplicateCount)
cts_clean <- cts_complete[!duplicated(cts_complete$ncbi_gene_id),]


# first column to rownames
cts_clean_rownames <- data.frame(cts_clean[,-1], row.names = cts_clean[,1])
```

## Gene ontology and pathway analysis

For this, column names must be in Entrez ID (=NCBI gene ID) format.

```{r}

##Enrichment in group 2 

# Gene ontoloy enrichment
#qlf <- glmQLFTest(et, coef=2)
go <- "goana"(et, species="Hs")
topGO(go, sort="up")

write.csv(as.data.frame(topGO(go, sort="up", n=20)), file= paste0(out_dir, "GOenrichment_up.csv"))

go <- "goana"(et, species="Hs")
topGO(go, sort="down")

write.csv(as.data.frame(topGO(go, sort="down", n=20)), file= paste0(out_dir, "GOenrichment_down.csv"))


