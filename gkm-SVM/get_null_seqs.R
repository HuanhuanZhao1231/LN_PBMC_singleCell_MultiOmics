# run_gkmSVM.R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a cell type as an argument.")

library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(data.table)

i <- args[1]
cat("Processing cell type:", i, "\n")

out_dir <- paste0("/public/home/zhaohuanhuan/snATAC/B/ArchR/gkmSVM/", i)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

peak_file <- paste0("~/snATAC/B/ArchR/realCellsubset/Save-proj_ProjHeme5/PeakCalls/PeakSet/", i, "_peaks.narrowPeak")
if (!file.exists(peak_file)) {
  stop(paste("File not found:", peak_file))
}
data <- fread(peak_file)
data <- data[1:20000, .(V1, V2, V3)]
data <- data[data$V1 %in% paste0("chr",1:22),]
bed_file <- paste0(out_dir, "/formated.", i, ".bed")
fwrite(data, bed_file, sep = "\t", col.names = FALSE)

genNullSeqs(paste0(out_dir, "/formated.", i, ".bed"),
  genomeVersion = 'hg19',
  outputBedFN = paste0(out_dir, "/negSet_", i, ".bed"),
  outputPosFastaFN = paste0(out_dir, "/posSet_", i, ".fa"),
  outputNegFastaFN = paste0(out_dir, "/negSet_", i, ".fa"),
  xfold = 1,
  repeat_match_tol = 0.02,
  GC_match_tol = 0.02,
  length_match_tol = 0.02,
  batchsize = 5000,
  nMaxTrials = 20,
  genome = NULL
)
