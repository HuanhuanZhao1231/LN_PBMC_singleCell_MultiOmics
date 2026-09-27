##########################################################################################
# Peak2Gene Links Methods
##########################################################################################

#' Add Peak2GeneLinks to an ArchRProject
#' 
#' This function will add peak-to-gene links to a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param useMatrix The name of the matrix containing gene expression information to be used for determining peak-to-gene links. See `getAvailableMatrices(ArchRProj)`
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param cellsToUse A character vector of cellNames to compute peak-to-gene links on if desired to run on a subset of the total cells.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from this analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param addEmpiricalPval Add empirical p-values based on randomly correlating peaks and genes not on the same seqname.
#' @param addPermutedPval Add permuted p-values based on shuffle sample correlating peaks and genes. This approach was adapted from
#' Regner et al 2021 "A multi-omic single-cell landscape of human gynecologic malignancies".
#' @param nperm An integer representing the number of permutations to run for Regner et al 2021 approach.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' 
#' @examples
#'
#' # Get Test ArchR Project
#' proj <- getTestProject()
#'
#' # Add P2G Links
#' proj <- addPeak2GeneLinks(proj, k = 20)
#'
#' # Get P2G Links
#' p2g <- getPeak2GeneLinks(proj)
#'
#' @export
addPeak2GeneLinks_P <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  excludeChr = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  addEmpiricalPval = FALSE,
  addPermutedPval = FALSE,
  nperm = 100,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
  ){
.validInput <- ArchR:::.validInput
.startLogging <- ArchR:::.startLogging
.endLogging <- ArchR:::.endLogging
.logThis <- ArchR:::.logThis
.logDiffTime <- ArchR:::.logDiffTime
.getFeatureDF <- ArchR:::.getFeatureDF
.computeKNN <- ArchR:::.computeKNN
.getGroupMatrix <- ArchR:::.getGroupMatrix
.getQuantiles <- ArchR:::.getQuantiles

.safelapply <- ArchR:::.safelapply
.safeSaveRDS <- ArchR:::.safeSaveRDS
.suppressAll <- ArchR:::.suppressAll

determineOverlapCpp <- ArchR:::determineOverlapCpp
rowCorCpp <- ArchR:::rowCorCpp

.validInput(
    input = ArchRProj,
    name = "ArchRProj",
    valid = c("ArchRProj")
)

.validInput(
    input = nperm,
    name = "nperm",
    valid = c("integer")
)


tstart <- Sys.time()

.startLogging(
    logFile = logFile
)


.logThis(
    mget(names(formals()), sys.frame(sys.nframe())),
    "addPermPeak2GeneLinks Input-Parameters",
    logFile = logFile
)
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = predictionCutoff, name = "predictionCutoff", valid = c("numeric", "null"))
  .validInput(input = addEmpiricalPval, name = "addEmpiricalPval", valid = c("boolean"))
  .validInput(input = addPermutedPval, name = "addPermutedPval", valid = c("boolean"))
  .validInput(input = nperm, name = "nperm", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeak2GeneLinks Input-Parameters", logFile = logFile)

  .logDiffTime(main="Getting Available Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)

  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }

  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }

  ArrowFiles <- getArrowFiles(ArchRProj)

  tstart <- Sys.time()

  dfAll <- .safelapply(seq_along(ArrowFiles), function(x){
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e){
      if(getArchRVerbose()) message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(
      cellNames = cNx,
      predictionScore = pSx
    )
  }, threads = threads) %>% Reduce("rbind", .)

  .logDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
    sum(dfAll[,2] < predictionCutoff), 
    nrow(dfAll), 
    round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), t1=tstart, verbose=verbose, logFile=logFile)

  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  if(!is.null(excludeChr)){
    peakSet <- peakSet[BiocGenerics::which(paste0(seqnames(peakSet)) %bcni% excludeChr)]
  }
  .logThis(peakSet, "peakSet", logFile = logFile)

  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  if(!is.null(excludeChr)){
    geneSet <- geneSet[BiocGenerics::which(geneSet$seqnames %bcni% excludeChr),,drop=FALSE]
  }
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  .logThis(geneStart, "geneStart", logFile = logFile)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)

  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)

  #Group Matrix RNA
  .logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
#    excludeSeqnames = excludeChr,
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  rawMatRNA <- groupMatRNA
  .logThis(groupMatRNA, "groupMatRNA", logFile = logFile)

  #Group Matrix ATAC
  .logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
#    excludeSeqnames = excludeChr,
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  rawMatATAC <- groupMatATAC
  .logThis(groupMatATAC, "groupMatATAC", logFile = logFile)

  .logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)

  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }

  names(geneStart) <- NULL

  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    rowRanges = geneStart
  )
  S4Vectors::metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)

  names(peakSet) <- NULL

  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    rowRanges = peakSet
  )
  S4Vectors::metadata(seATAC)$KNNList <- knnObj
  .logThis(seATAC, "seATAC", logFile = logFile)

  rm(groupMatRNA, groupMatATAC)
  gc()

  #Overlaps
  .logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
  o <- DataFrame(
    findOverlaps(
      .suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )

  #Get Distance from Fixed point A B 
  o$distance <- GenomicRanges::distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")

  #Compute PVal Stats
  .logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "Pval", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "Pval", "FDR", "VarQATAC", "VarQRNA")  
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  S4Vectors::metadata(out)$peakSet <- peakSet
  S4Vectors::metadata(out)$geneSet <- geneStart

  #Null Correlations
  if(addEmpiricalPval){
    .logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
    nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
    out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }

  #Permuted Pval
  if(addPermutedPval){
    message("Performing Permuted P-values similar to Regner et al., 2021")
    #Permute
    p <- o
    o$PermPval <- 0
    for(i in seq_len(nperm)){
      message("Running Permutation ", i, " of ", nperm)
      idx <- sample(ncol(seATAC))
      p$Correlation <- rowCorCpp(as.integer(p$A), as.integer(p$B), assay(seATAC)[,idx,drop=FALSE], assay(seRNA))
      p$TStat <- (p$Correlation / sqrt((1-p$Correlation^2)/(ncol(seATAC)-2))) #T-statistic P-value
      p$Pval <- 2*pt(-abs(p$TStat), ncol(seATAC) - 2)
      cdf <- ecdf(p$Pval)
      o$PermPval <- o$PermPval + p$Pval
    }
    o$PermPval <- (o$PermPval / nperm) #Average
    o$PermFDR <- pmin(ecdf(o$PermPval)(o$Pval) / ecdf(o$Pval)(o$Pval), 1)
  }

  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN_Pval.rds")
  .safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN_Pval.rds")
  .safeSaveRDS(seRNA, outRNA, compress = FALSE)
  S4Vectors::metadata(out)$seATAC <- outATAC
  S4Vectors::metadata(out)$seRNA <- outRNA

  S4Vectors::metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out

  .logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  .endLogging(logFile = logFile)

  ArchRProj

}

.getNullCorrelations <- function(seA, seB, o, n){

  o$seq <- seqnames(seA)[o$A]

  nullCor <- lapply(seq_along(unique(o$seq)), function(i){

    #Get chr from olist
    chri <- unique(o$seq)[i]
    #message(chri, " ", appendLF = FALSE)

    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)

    return(out)

  }) %>% SimpleList
  #message("")

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

}
