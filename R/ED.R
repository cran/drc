"ED" <-
function(obj, percVec, od = FALSE, ...)
{
#    typeStr <- "ED"
#
#
#    ## Finding super class
#    EDinfo <- partIn(class(obj), typeStr)
#
#    EDstr <- EDinfo[[1]]
#    EDfct <- EDinfo[[2]]
#
#    inOut <- iofct(typeStr, EDstr)
#    in1fct <- inOut[[1]]
#    in2fct <- inOut[[2]]
#    outfct <- inOut[[3]]


#    EDlist <- obj[[11]][[8]]

    if (inherits(obj, "bindrc"))
    {
        ED2(obj, percVec)
    } else {

    ## Checking 'percVec' vector ... should be numbers between 0 and 100
    if (any(percVec<=0 | percVec>=100)) {stop("Percentages outside the interval ]0, 100[ not allowed")}
    lenPV <- length(percVec)
    

    ## Retrieving relevant quantities
    EDlist <- obj$fct$"edfct"    
    assayNo <- obj$data[, 3]  # obj[[9]][,3]
    numAss <- length(unique(assayNo))

    parmMat <- obj$"parmMat"  # [[10]]
    sumObj <- summary(obj, od = od)
    resVar <- sumObj$"resVar"  # [[1]]
    varMat<-obj$"transformation"%*%sumObj$"varMat"%*%t(obj$"transformation")        
#    varMat <- obj[[12]]%*%sumObj[[2]]%*%t(obj[[12]])
#    varMat <- sumObj[[2]]
#    parm <- c((sumObj[[3]])[,1])
    parm <- c((sumObj$"estimates")[,1])
#    strParm <- (unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2] 
#    strParm <- unique(obj[[9]][, ncol(obj[[9]]) - 1])  # second last column contains original curve levels
    strParm <- unique(obj$"data"[, ncol(obj$"data") - 1])  # second last column contains original curve levels
#    print(strParm)
#    strParm <- strParm[apply(parmMat, 2, function(x){!any(is.na(x))})]

#    ncPM <- ncol(parmMat)
#    naVec <- rep(NA, ncPM)
#    for (i in 1:ncPM)
#    {
#        if (any(is.na(parmMat[,i]))) {naVec[i] <- i}
#    }


    ## Creating an index matrix
#    asVec1 <- c(t(parmMat))
#    notNA <- !is.na(asVec1)
#    asVec2 <- asVec1[notNA]
#    asVec1[notNA] <- match(asVec2, unique(asVec2))
#    indexMat <- matrix(asVec1, nrow(parmMat), ncol(parmMat), byrow=TRUE)

#    lenNV <- length(naVec[!is.na(naVec)])
#    if (lenNV>0)
#    {
#        parmMat <- parmMat[,-naVec[!is.na(naVec)]]
#        indexMat <- indexMat[,-naVec[!is.na(naVec)]]
#    }
    
    ncPM <- ncol(parmMat)  # obj[[10]])
    nrPM <- nrow(parmMat)  # obj[[10]])
#    stop()
    
#    options(warn=-1)  # to avoid warnings when filling in matrix with elements in excess in the vector 
#    indexMat <- t(matrix(NA, nrPM, ncPM))
#    indexMat[!is.na(t(obj[[10]]))] <- 1:(nrPM*ncPM)
#    indexMat[!is.na(t(parmMat))] <- 1:(nrPM*ncPM)
#    indexMat <- t(indexMat)
#    options(warn=0)
#    print(indexMat)
#    stop()
    
    naVec <- rep(FALSE, ncPM)
    for (i in 1:ncPM)
    {
        naVec[i] <- any(is.na(parmMat[, i]))
    }
#    indexMat <- indexMat[, !naVec, drop=FALSE]
    parmMat <- parmMat[, !naVec, drop=FALSE] 
    strParm <- strParm[!naVec]
#    print(indexMat)

    ncPM2 <- ncol(parmMat)  # obj[[10]])
    nrPM2 <- nrow(parmMat)  # obj[[10]])
    indexMat <- matrix(1:(nrPM2*ncPM2), nrPM2, ncPM2, byrow = TRUE)   


    ## Finding out which parameter occurs most times; this determines the number of ED values
#    maxIndex <- 0
#    maxParm <- 0
#    
#    print(EDlist(parmMat[,1], 50, upper))
#    indexVec <- (1:nrow(parmMat))[(EDlist(parmMat[,1], 50, upper))[[2]]<1e-10]
#    print(indexVec)
#    for (i in indexVec)
#    {
#        PM <- parmMat[i,]
#        lenPM <- length(unique(PM))
#        if (lenPM > maxParm) {maxIndex <- match(unique(PM),PM); maxParm <- i}
#    }
#
#    nCol <- ncol(parmMat)
#    indexVec <- 1:nCol
#    for (i in 1:nCol)
#    {
#        if (any(is.na(parmMat[,i]))) {indexVec[i] <- NA}
#    }
#    indexVec <- indexVec[!is.na(indexVec)]
#    print(indexVec)
    indexVec <- 1:ncol(indexMat)


    ## Calculating ED values    
    lenEB <- length(indexVec)    
    dimNames <- rep("", lenEB*lenPV)
    EDmat <- matrix(0, lenEB*lenPV, 2)

#print(dim(EDmat))
#print(parmMat)
#print(indexVec)
#print(dimNames)
#print(indexMat)    
    
    rowIndex <- 1
#    for (i in maxIndex)
    for (i in indexVec)
    {
        parmInd <- indexMat[, i]
        varCov <- varMat[parmInd, parmInd]
        parmChosen <- parmMat[, i]

        for (j in 1:lenPV)
        {
            EDeval <- EDlist(parmChosen, percVec[j], ...)
            EDmat[rowIndex, 1] <- EDeval[[1]]
            dfEval <- EDeval[[2]]
#            print(dfEval)
        
            EDmat[rowIndex, 2] <- sqrt(dfEval%*%varCov%*%dfEval)

            dimNames[rowIndex] <- paste(strParm[i], ":", percVec[j], sep="")
            rowIndex <- rowIndex + 1
        }
    }
    dimnames(EDmat) <- list(dimNames, c("Estimate", "Std. Error"))

    ## Displaying the ED values
    cat("\n")
    printCoefmat(EDmat)
    invisible(EDmat)    
    }
}
