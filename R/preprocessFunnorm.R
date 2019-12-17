######################################################
## Functional normalization of the 450k array
## Jean-Philippe Fortin
## Sept 24 2013
#####################################################

## Return a normalized beta matrix
##
## need stuff from gmodels.  Where?
preprocessFunnorm <- function(rgSet, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE, verbose = TRUE) {

    .isMatrixBackedOrStop(rgSet, "preprocessFunnorm")

    .isRGOrStop(rgSet)
    rgSet <- updateObject(rgSet) ## FIXM: might not KDH: technically, this should not be needed, but might be nice

    # Background correction and dye bias normalization:
    if (bgCorr){
        if(verbose && dyeCorr) {
            message("[preprocessFunnorm] Background and dye bias correction with noob")
        } else {
            message("[preprocessFunnorm] Background correction with noob")
        }
        invisible(gc())
        gmSet <- preprocessNoob(rgSet, dyeCorr = dyeCorr, verbose=verbose)
        if(verbose) message("[preprocessFunnorm] Mapping to genome")
        invisible(gc())
        gmSet <- mapToGenome(gmSet)
    } else {
        if(verbose) message("[preprocessFunnorm] Mapping to genome")
        invisible(gc())
        gmSet <- mapToGenome(rgSet)
    }

    subverbose <- max(as.integer(verbose) - 1L, 0)

    if(verbose) message("[preprocessFunnorm] Quantile extraction")
    invisible(gc())
    extractedData <- .extractFromRGSet450k(rgSet)
    rm(rgSet)
    invisible(gc())

    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    if(verbose) message("[preprocessFunnorm] Normalization...")
    if(keepCN) {
        if(verbose) message("[preprocessFunnorm] Normalization-getCN")
        CN <- getCN(gmSet)
    }
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k")
    gmSet <- .normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                   sex = sex, nPCs = nPCs, verbose = subverbose)
    preprocessMethod <- c(preprocessMethod(gmSet),
                          mu.norm = sprintf("Funnorm, nPCs=%s", nPCs))
    if(ratioConvert) {
        invisible(gc())
        if(verbose) message("[preprocessFunnorm] Normalization-ratioConvert")
        grSet <- ratioConvert(gmSet, type = "Illumina", keepCN = keepCN)
        if(keepCN) {
            assay(grSet, "CN") <- CN
        }
        grSet@preprocessMethod <- preprocessMethod
        return(grSet)
    } else {
        gmSet@preprocessMethod <- preprocessMethod
        return(gmSet)
    }
    invisible(gc())
 }

 .getFunnormIndices <- function(object) {
     .isGenomicOrStop(object)
     probeType <- getProbeType(object, withColor = TRUE)
     autosomal <- (seqnames(object) %in% paste0("chr", 1:22))
     indices <- list(IGrn = which(probeType == "IGrn" & autosomal),
                     IRed = which(probeType == "IRed" & autosomal),
                     II = which(probeType == "II" & autosomal),
                     X = which(seqnames(object) == "chrX"),
                     Y = which(seqnames(object) == "chrY"))
     indices
 }

 .normalizeFunnorm450k <- function(object, extractedData, nPCs, sex, verbose = TRUE) {
     normalizeQuantiles <- function(matrix, indices, sex = NULL, verbose=TRUE) {
         matrix <- matrix[indices,,drop=FALSE]
         ## uses probs, model.matrix, nPCS, through scoping)
         if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-sex")
         oldQuantiles <- t(colQuantiles(matrix, probs = probs))
         if(is.null(sex)) {
             newQuantiles <- .returnFit(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs)
         } else {
             newQuantiles <- .returnFitBySex(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs, sex = sex)
         }
         if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.normalizeMatrix")
         .normalizeMatrix(matrix, newQuantiles, verbose=verbose)
     }
    
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl")
    indicesList <- .getFunnormIndices(object)
    model.matrix <- .buildControlMatrix450k(extractedData, verbose=verbose)
    probs <- seq(from = 0, to = 1, length.out = 500)
     
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-getMeth/getUnmeth")
    Meth <- getMeth(object)
    Unmeth <- getUnmeth(object)
    invisible(gc())
    if (nPCs > 0){
        if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-nPCs")
        for (type in c("IGrn", "IRed", "II")) {
            invisible(gc())
            indices <- indicesList[[type]]
            if(length(indices) > 0) {
                if(verbose) message(sprintf("[normalizeFunnorm450k] Normalization of the %s probes", type))
                Unmeth[indices,] <- normalizeQuantiles(Unmeth, indices = indices, sex = NULL, verbose=verbose)
                Meth[indices,] <- normalizeQuantiles(Meth, indices = indices, sex = NULL, verbose=verbose)
            }
        }

        invisible(gc())
        if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-X")
        indices <- indicesList[["X"]]
        if(length(indices) > 0) {
            if(verbose) message("[normalizeFunnorm450k] Normalization of the X-chromosome")
            Unmeth[indices,] <- normalizeQuantiles(Unmeth, indices = indices, sex = sex, verbose=verbose)
            Meth[indices,] <- normalizeQuantiles(Meth, indices = indices, sex = sex, verbose=verbose)
        }
    }

    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-Y")
    indices <- indicesList[["Y"]]
    if(length(indices) > 0) {
        if(verbose) message("[normalizeFunnorm450k] Normalization of the Y-chromosome")
        sex <- as.character(sex)
        levels <- unique(sex)
        nSexes <- length(levels)
        if (nSexes == 2) {
            level1 <- levels[1]
            level2 <- levels[2]
        }
        if (nSexes == 2) {
            if (sum(sex == level1)>1) {
                invisible(gc())
                Meth[indices, sex==level1]   <- preprocessCore::normalize.quantiles(Meth[indices, sex == level1, drop=FALSE])
                Unmeth[indices, sex==level1] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level1,drop=FALSE])
            }
            if (sum(sex == level2)>1) {
                invisible(gc())
                Meth[indices, sex==level2]   <- preprocessCore::normalize.quantiles(Meth[indices, sex == level2,drop=FALSE])
                Unmeth[indices, sex==level2] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level2,drop=FALSE])
            }
        } else {
            invisible(gc())
            Meth[indices,] <- preprocessCore::normalize.quantiles(Meth[indices,])
            Unmeth[indices,] <- preprocessCore::normalize.quantiles(Unmeth[indices,])
        }
    }
    assay(object, "Meth") <- Meth
    assay(object, "Unmeth") <- Unmeth
    return(object)
}



### To extract quantiles and control probes from rgSet
.extractFromRGSet450k <- function(rgSet) {
    invisible(gc())
    rgSet <- updateObject(rgSet)
    controlType <- c("BISULFITE CONVERSION I",
                     "BISULFITE CONVERSION II",
                     "EXTENSION",
                     "HYBRIDIZATION",
                     "NEGATIVE",
                     "NON-POLYMORPHIC",
                     "NORM_A",
                     "NORM_C",
                     "NORM_G",
                     "NORM_T",
                     "SPECIFICITY I",
                     "SPECIFICITY II",
                     "TARGET REMOVAL",
                     "STAINING")

    array <- annotation(rgSet)[["array"]]
    ## controlAddr <- getControlAddress(rgSet, controlType = controlType, asList = TRUE)
    ctrls <- getProbeInfo(rgSet, type = "Control")
    if(!all(controlType %in% ctrls$Type))
        stop("The `rgSet` does not contain all necessary control probes")
    invisible(gc())
    ctrlsList <- split(ctrls, ctrls$Type)[controlType]
    redControls <- getRed(rgSet)[ctrls$Address,,drop=FALSE]
    redControls <- lapply(ctrlsList, function(ctl) redControls[ctl$Address,,drop=FALSE])
    greenControls <- getGreen(rgSet)[ctrls$Address,,drop=FALSE]
    greenControls <- lapply(ctrlsList, function(ctl) greenControls[ctl$Address,,drop=FALSE])

    ## Extraction of the undefined negative control probes
    oobRaw <- getOOB(rgSet)
    probs <- c(0.01, 0.50, 0.99)
    greenOOB <- t(colQuantiles(oobRaw$Grn, na.rm = TRUE, probs = probs))
    redOOB   <- t(colQuantiles(oobRaw$Red, na.rm=TRUE,  probs = probs))
    oob      <- list(greenOOB = greenOOB, redOOB = redOOB)
    invisible(gc())

    return(list(
        greenControls = greenControls,
        redControls = redControls,
        oob = oob, ctrlsList = ctrlsList,
        array = array))
}


## Extraction of the Control matrix
.buildControlMatrix450k <- function(extractedData, verbose=TRUE) {
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-setVars")
    getCtrlsAddr <- function(exType, index) {
        ctrls <- ctrlsList[[index]]
        addr <- ctrls$Address
        names(addr) <- ctrls$ExtendedType
        na.omit(addr[exType])
    }

    array <- extractedData$array
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)
    ctrlsList <- extractedData$ctrlsList

    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-BC2")
    ## Bisulfite conversion extraction for probe type II:
    index <- match("BISULFITE CONVERSION II", controlNames)
    redControls.current <- redControls[[ index ]]
    bisulfite2 <- colMeans2(redControls.current, na.rm = TRUE)

    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-BC1")
    ## Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I", controlNames)
    if (array=="IlluminaHumanMethylation450k"){
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3), index = index)
    } else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c("-", "-"), 1:2), index = index)
    }
    greenControls.current <- greenControls[[ index ]][addr,,drop=FALSE]
    if (array=="IlluminaHumanMethylation450k"){
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 4:6), index = index)
    } else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 3:5), index = index)
    }
    redControls.current <- redControls[[ index ]][addr,, drop=FALSE]
    if (nrow(redControls.current)==nrow(greenControls.current)){
        bisulfite1 <- colMeans2(redControls.current + greenControls.current, na.rm = TRUE)
    } else {
        bisulfite1 <- colMeans2(redControls.current, na.rm=TRUE) + colMeans2(greenControls.current, na.rm = TRUE)
    }


    ## Staining
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Staining")
    invisible(gc())
    index <- match("STAINING", controlNames)
    addr <- getCtrlsAddr(exType = "Biotin (High)", index = index)
    stain.green <- t(greenControls[[ index ]][addr,,drop=FALSE])
    addr <- getCtrlsAddr(exType = "DNP (High)", index = index)
    stain.red <- t(redControls[[ index ]][addr,, drop=FALSE ])

    ## Extension
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Extension")
    invisible(gc())
    index <-    match("EXTENSION", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("A", "T")), index = index)
    extension.red <- t(redControls[[index]][addr,,drop=FALSE])
    colnames(extension.red) <- paste0("extRed", 1:ncol(extension.red))
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("C", "G")), index = index)
    extension.green <- t(greenControls[[index]][addr,,drop=FALSE])
    colnames(extension.green) <- paste0("extGrn", 1:ncol(extension.green))

    ## Hybridization should be monitored only in the green channel
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Hybridisation")
    index <- match("HYBRIDIZATION", controlNames)
    hybe <- t(greenControls[[index]])
    colnames(hybe) <- paste0("hybe", 1:ncol(hybe))

    ## Target removal should be low compared to hybridization probes
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Targetrem")
    index <- match("TARGET REMOVAL", controlNames)
    targetrem <- t(greenControls[[index]])
    colnames(targetrem) <- paste0("targetrem", 1:ncol(targetrem))

    ## Non-polymorphic probes
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Non-polymorphic-probes")
    invisible(gc())
    index <- match("NON-POLYMORPHIC", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("A", "T")), index = index)
    nonpoly.red <- t(redControls[[index]][addr, ,drop=FALSE])
    colnames(nonpoly.red) <- paste0("nonpolyRed", 1:ncol(nonpoly.red))
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("C", "G")), index = index)
    nonpoly.green <- t(greenControls[[index]][addr, ,drop=FALSE])
    colnames(nonpoly.green) <- paste0("nonpolyGrn", 1:ncol(nonpoly.green))

    ## Specificity II
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Specificity2")
    invisible(gc())
    index <- match("SPECIFICITY II", controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    spec2.green <- t(greenControls.current)
    colnames(spec2.green) <- paste0("spec2Grn", 1:ncol(spec2.green))
    spec2.red <- t(redControls.current)
    colnames(spec2.red) <- paste0("spec2Red", 1:ncol(spec2.red))
    spec2.ratio <- colMeans2(greenControls.current, na.rm = TRUE) /
        colMeans2(redControls.current, na.rm = TRUE)

    ## Specificity I
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Specificity1")
    invisible(gc())
    index <- match("SPECIFICITY I", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 1:3), index = index)
    greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
    redControls.current <- redControls[[index]][addr,,drop=FALSE]
    spec1.green <- t(greenControls.current)
    colnames(spec1.green) <- paste0("spec1Grn", 1:ncol(spec1.green))
    spec1.ratio1 <- colMeans2(redControls.current, na.rm = TRUE) /
        colMeans2(greenControls.current, na.rm = TRUE)

    index <- match("SPECIFICITY I", controlNames) # Added that line
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 4:6), index = index)
    greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
    redControls.current <- redControls[[index]][addr,,drop=FALSE]
    spec1.red <- t(redControls.current)
    colnames(spec1.red) <- paste0("spec1Red", 1:ncol(spec1.red))
    spec1.ratio2 <- colMeans2(greenControls.current, na.rm = TRUE) /
        colMeans2(redControls.current, na.rm = TRUE)
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2) / 2

    ## Normalization probes:
    invisible(gc())
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Normalisation_Probes")
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans2(greenControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans2(greenControls[[index]], na.rm = TRUE)

    dyebias <- (normC + normG)/(normA + normT)

    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    oob.ratio <- oobG[2,]/oobR[2,]
    oobG <- t(oobG)
    colnames(oobG) <- paste0("oob", c(1,50,99))
    
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-cbind")
    invisible(gc())
    model.matrix <- cbind(
        bisulfite1, bisulfite2, extension.green, extension.red, hybe,
        stain.green, stain.red, nonpoly.green, nonpoly.red,
        targetrem, spec1.green, spec1.red, spec2.green, spec2.red, spec1.ratio1,
        spec1.ratio, spec2.ratio, spec1.ratio2, normA, normC, normT, normG, dyebias,
        oobG, oob.ratio)


    ## Imputation
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Imputation")
    invisible(gc())
    for (colindex in 1:ncol(model.matrix)) {
        if(any(is.na(model.matrix[,colindex]))) {
            column <- model.matrix[,colindex]
            column[is.na(column)]    <- mean(column, na.rm = TRUE)
            model.matrix[ , colindex] <- column
        }
    }

    ## Scaling
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Scaling")
    invisible(gc())
    model.matrix <- scale(model.matrix)

    ## Fixing outliers
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Outliers")
    invisible(gc())
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3

    ## Rescaling
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.buildControl-Rescaling")
    invisible(gc())
    model.matrix <- scale(model.matrix)

    invisible(gc())
    return(model.matrix)
}


### Return the normalized quantile functions
.returnFit <- function(controlMatrix, quantiles, nPCs) {
    invisible(gc())
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    ## Fixing potential problems with extreme quantiles
    quantiles[1,] <- 0
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles) - 1,] + 1000
    meanFunction <- rowMeans2(quantiles)
    res <- quantiles - meanFunction
    controlPCs <- prcomp(controlMatrix)$x[,1:nPCs,drop=FALSE]
    design <- model.matrix(~controlPCs)
    fits <- lm.fit(x = design, y = t(res))
    newQuantiles <- meanFunction + t(fits$residuals)
    newQuantiles <- .regularizeQuantiles(newQuantiles)
    invisible(gc())
    return(newQuantiles)
}

.returnFitBySex <- function(controlMatrix, quantiles, nPCs, sex) {
    invisible(gc())
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    sex    <- as.character(sex)
    levels <- unique(sex)
    nSexes <- length(levels)
    if (nSexes == 2) {
        sex1 <- sum(sex == levels[1])
        sex2 <- sum(sex == levels[2])

    } else {
        sex1 <- sum(sex == levels[1])
        sex2 <- 0
    }

    ## When normalization should not be performed by sex separately:
    if ((sex1 <= 10) | (sex2 <= 10)) {
        newQuantiles <- .returnFit(controlMatrix = controlMatrix,
                                  quantiles = quantiles,
                                  nPCs = nPCs)
    } else {
        quantiles1 <- quantiles[, sex == levels[1]]
        controlMatrix1 <- controlMatrix[sex == levels[1], ]

        newQuantiles1 <- .returnFit(controlMatrix = controlMatrix1,
                                   quantiles = quantiles1,
                                   nPCs = nPCs)

        quantiles2 <- quantiles[, sex == levels[2]]
        controlMatrix2 <- controlMatrix[sex == levels[2], ]

        newQuantiles2 <- .returnFit(controlMatrix = controlMatrix2,
                                   quantiles = quantiles2,
                                   nPCs = nPCs)

        newQuantiles <- quantiles
        newQuantiles[, sex == levels[1]] <- newQuantiles1
        newQuantiles[, sex == levels[2]] <- newQuantiles2
    }

    invisible(gc())
    return(newQuantiles)
}

### Normalize a matrix of intensities
.normalizeMatrix <- function(intMatrix, newQuantiles, verbose=TRUE) {
    invisible(gc())
    ## normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
    n <- nrow(newQuantiles)
    if(verbose) message("[preprocessFunnorm] Normalization-.normalizeFunnorm450k-.normalizeMatrix-sapply")
    normMatrix <- sapply(1:ncol(intMatrix), function(i) {
        crtColumn <- intMatrix[ , i]
        crtColumn.reduced <- crtColumn[!is.na(crtColumn)]
        ## Generation of the corrected intensities:
        target <- sapply(1:(n-1), function(j) {
            start <- newQuantiles[j,i]
            end <- newQuantiles[j+1,i]
            if (!isTRUE(all.equal(start,end))){
                sequence <- seq(start, end,( end-start)/n)[-(n+1)]
            } else {
                sequence <- rep(start, n)
            }
            return(sequence)
        })
        target <- as.vector(target)
        result <- preprocessCore::normalize.quantiles.use.target(matrix(crtColumn.reduced), target)
        return(result)
    })
    invisible(gc())
    return(normMatrix)
}

# To ensure a monotonically increasing and non-negative quantile function
# Necessary for pathological cases
.regularizeQuantiles <- function(x){
    x[x<0] <- 0
    colCummaxs(x)
}
