addPermPeak2GeneLinks <- function(
    ArchRProj = NULL,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 1:30,
    scaleDims = NULL,
    corCutOff = 0.5,
    cellsToUse = NULL,
    k = 100,
    knnIteration = 500,
    overlapCutoff = 0.8,
    maxDist = 250000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    predictionCutoff = 0.4,
    addEmpiricalPval = FALSE,
    addPermutedPval = FALSE,
    nperm = 1000,
    seed = 1,
    threads = max(floor(getArchRThreads()/2),1),
    verbose = TRUE,
    logFile = createLogFile("addPermPeak2GeneLinks")
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


############################################################
## Check matrices
############################################################

.logDiffTime(
    main="Getting Available Matrices",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)


AvailableMatrices <- getAvailableMatrices(ArchRProj)


if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
}


if(useMatrix %ni% AvailableMatrices){
    stop(
        paste0(
            useMatrix,
            " not in AvailableMatrices"
        )
    )
}



############################################################
## Prediction score filtering
############################################################


ArrowFiles <- getArrowFiles(ArchRProj)


dfAll <- .safelapply(
    seq_along(ArrowFiles),
    function(x){

        cNx <- paste0(
            names(ArrowFiles)[x],
            "#",
            h5read(
                ArrowFiles[x],
                paste0(
                    useMatrix,
                    "/Info/CellNames"
                )
            )
        )


        pSx <- tryCatch(
            {
                h5read(
                    ArrowFiles[x],
                    paste0(
                        useMatrix,
                        "/Info/predictionScore"
                    )
                )
            },
            error=function(e){
                rep(
                    9999999,
                    length(cNx)
                )
            }
        )


        DataFrame(
            cellNames=cNx,
            predictionScore=pSx
        )

    },
    threads=threads
) %>% Reduce("rbind", .)



keep <- sum(
    dfAll[,2] >= predictionCutoff
) /
nrow(dfAll)



dfAll <- dfAll[
    which(
        dfAll[,2] > predictionCutoff
    ),
]



############################################################
## Peak and gene annotation
############################################################


set.seed(seed)


peakSet <- getPeakSet(ArchRProj)



geneSet <- .getFeatureDF(
    ArrowFiles,
    useMatrix,
    threads=threads
)



geneStart <- GRanges(
    geneSet$seqnames,
    IRanges(
        geneSet$start,
        width=1
    ),
    name=geneSet$name,
    idx=geneSet$idx
)



############################################################
## Reduced dimension and KNN groups
############################################################


rD <- getReducedDims(
    ArchRProj,
    reducedDims=reducedDims,
    corCutOff=corCutOff,
    dimsToUse=dimsToUse
)


if(!is.null(cellsToUse)){
    rD <- rD[
        cellsToUse,
        ,
        drop=FALSE
    ]
}



idx <- sample(
    seq_len(nrow(rD)),
    knnIteration,
    replace=!nrow(rD)>=knnIteration
)



.logDiffTime(
    main="Computing KNN",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)



knnObj <- .computeKNN(
    data=rD,
    query=rD[idx,],
    k=k
)



keepKnn <- determineOverlapCpp(
    knnObj,
    floor(overlapCutoff*k)
)



knnObj <- knnObj[
    keepKnn==0,
]



knnObj <- lapply(
    seq_len(nrow(knnObj)),
    function(x){
        rownames(rD)[knnObj[x,]]
    }
) %>% SimpleList



############################################################
## Group matrices
############################################################


geneDF <- mcols(geneStart)
peakDF <- mcols(peakSet)

geneDF$seqnames <- seqnames(geneStart)
peakDF$seqnames <- seqnames(peakSet)



.logDiffTime(
    main="Getting Group RNA Matrix",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)



groupMatRNA <- .getGroupMatrix(
    ArrowFiles=getArrowFiles(ArchRProj),
    featureDF=geneDF,
    groupList=knnObj,
    useMatrix=useMatrix,
    threads=threads,
    verbose=FALSE
)



rawMatRNA <- groupMatRNA



.logDiffTime(
    main="Getting Group ATAC Matrix",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)



groupMatATAC <- .getGroupMatrix(
    ArrowFiles=getArrowFiles(ArchRProj),
    featureDF=peakDF,
    groupList=knnObj,
    useMatrix="PeakMatrix",
    threads=threads,
    verbose=FALSE
)



rawMatATAC <- groupMatATAC



############################################################
## Normalize
############################################################


groupMatRNA <- 
    t(
        t(groupMatRNA) /
        colSums(groupMatRNA)
    ) *
    scaleTo



groupMatATAC <-
    t(
        t(groupMatATAC) /
        colSums(groupMatATAC)
    ) *
    scaleTo



if(log2Norm){

    groupMatRNA <- log2(groupMatRNA+1)

    groupMatATAC <- log2(groupMatATAC+1)

}



############################################################
## Create SummarizedExperiment
############################################################


names(geneStart)<-NULL


seRNA <- SummarizedExperiment(
    assays=
        SimpleList(
            RNA=groupMatRNA,
            RawRNA=rawMatRNA
        ),
    rowRanges=geneStart
)


metadata(seRNA)$KNNList <- knnObj



names(peakSet)<-NULL


seATAC <- SummarizedExperiment(
    assays=
        SimpleList(
            ATAC=groupMatATAC,
            RawATAC=rawMatATAC
        ),
    rowRanges=peakSet
)


metadata(seATAC)$KNNList <- knnObj
############################################################
## Find Peak-Gene pairs
############################################################


.logDiffTime(
    main="Finding Peak Gene Pairings",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)


o <- DataFrame(
    findOverlaps(
        .suppressAll(
            resize(
                seRNA,
                2 * maxDist + 1,
                "center"
            )
        ),
        resize(
            rowRanges(seATAC),
            1,
            "center"
        ),
        ignore.strand=TRUE
    )
)



o$distance <- distance(
    rowRanges(seRNA)[o[,1]],
    rowRanges(seATAC)[o[,2]]
)


colnames(o) <- c(
    "B",
    "A",
    "distance"
)



############################################################
## Observed correlation
############################################################


.logDiffTime(
    main="Computing Observed Correlations",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)



o$Correlation <- rowCorCpp(
    as.integer(o$A),
    as.integer(o$B),
    assay(seATAC),
    assay(seRNA)
)



o$VarAssayA <- .getQuantiles(
    matrixStats::rowVars(
        assay(seATAC)
    )
)[o$A]


o$VarAssayB <- .getQuantiles(
    matrixStats::rowVars(
        assay(seRNA)
    )
)[o$B]



############################################################
## Original ArchR FDR
############################################################


o$TStat <- (
    o$Correlation /
        sqrt(
            pmax(
                1 - o$Correlation^2,
                1e-17,
                na.rm=TRUE
            ) /
                (ncol(seATAC)-2)
        )
)



o$Pval <- 2 *
    pt(
        -abs(o$TStat),
        ncol(seATAC)-2
    )


o$FDR <- p.adjust(
    o$Pval,
    method="fdr"
)



############################################################
## Permutation-based null model
############################################################


if(addPermutedPval){

    .logDiffTime(
        main=paste0(
            "Computing Permutation Null Correlations (",
            nperm,
            " permutations)"
        ),
        t1=tstart,
        verbose=verbose,
        logFile=logFile
    )


    obsCor <- o$Correlation


    nullExtreme <- numeric(
        length(obsCor)
    )


    seATAC.mat <- assay(seATAC)



    for(i in seq_len(nperm)){


        if(verbose){

            message(
                "Permutation ",
                i,
                "/",
                nperm
            )

        }



        ## Randomly shuffle ATAC aggregate labels
        permATAC <- seATAC.mat[
            ,
            sample(
                ncol(seATAC.mat)
            ),
            drop=FALSE
        ]



        nullCor <- rowCorCpp(
            as.integer(o$A),
            as.integer(o$B),
            permATAC,
            assay(seRNA)
        )



        nullExtreme <- nullExtreme +
            (
                nullCor >= obsCor
            )


    }



    ## empirical permutation p value

    o$PermPval <- (
        nullExtreme + 1
    ) /
        (
            nperm + 1
        )



    o$PermFDR <- p.adjust(
        o$PermPval,
        method="BH"
    )


}



############################################################
## Output formatting
############################################################


keepCols <- c(
    "A",
    "B",
    "Correlation",
    "FDR",
    "VarAssayA",
    "VarAssayB"
)



if(addPermutedPval){

    keepCols <- c(
        keepCols,
        "PermPval",
        "PermFDR"
    )

}



out <- o[,keepCols]



colnames(out)[1:6] <- c(
    "idxATAC",
    "idxRNA",
    "Correlation",
    "FDR",
    "VarQATAC",
    "VarQRNA"
)



############################################################
## Save metadata
############################################################


mcols(peakSet) <- NULL
names(peakSet) <- NULL


metadata(out)$peakSet <- peakSet
metadata(out)$geneSet <- geneStart



dir.create(
    file.path(
        getOutputDirectory(ArchRProj),
        "Peak2GeneLinks"
    ),
    showWarnings=FALSE
)



outATAC <- file.path(
    getOutputDirectory(ArchRProj),
    "Peak2GeneLinks",
    "seATAC-Group-KNN.rds"
)



.safeSaveRDS(
    seATAC,
    outATAC,
    compress=FALSE
)



outRNA <- file.path(
    getOutputDirectory(ArchRProj),
    "Peak2GeneLinks",
    "seRNA-Group-KNN.rds"
)



.safeSaveRDS(
    seRNA,
    outRNA,
    compress=FALSE
)



metadata(out)$seATAC <- outATAC
metadata(out)$seRNA <- outRNA



############################################################
## Store Peak2GeneLinks
############################################################


metadata(
    ArchRProj@peakSet
)$Peak2GeneLinks <- out



.logDiffTime(
    main="Completed Permutation Peak2Gene Correlations!",
    t1=tstart,
    verbose=verbose,
    logFile=logFile
)



.endLogging(
    logFile=logFile
)



return(ArchRProj)

}
