"SI" <-
function(obj, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, ...)
{
#    typeStr <- "SI"
#
#
#    ## Finding super class
#    SIinfo <- partIn(class(obj), typeStr)
#
#    SIstr <- SIinfo[[1]]
#    SIfct <- SIinfo[[2]]
#
#    inOut <- iofct(typeStr, SIstr)
#    in1fct <- inOut[[1]]
#    in2fct <- inOut[[2]]
#    outfct <- inOut[[3]]


#    SIlist <- obj[[11]][[9]]
    SIlist <- obj$fct$"sifct"
    if (is.null(SIlist)) {stop("SI values cannot be calculated")}

    ## Checking contain of percVec vector ... should be numbers between 0 and 100
    if (any(percVec<=0 | percVec>=100)) {stop("Percentages outside the interval [0, 100] not allowed")}

    if (missing(compMatch)) {matchNames <- FALSE} else {matchNames <- TRUE}


    ## Retrieving relevant quantities
    assayNo <- obj$"data"[, 3]  # [[9]][,3]
    numAss <- length(unique(assayNo))

    parmMat <- obj$"parmMat"  # [[10]]
    sumObj<-summary(obj, od = od)
    varMat<-obj$"transformation"%*%sumObj$"varMat"%*%t(obj$"transformation")    
#    varMat<-obj[[12]]%*%sumObj[[2]]%*%t(obj[[12]])
#    parm<-c((sumObj[[3]])[,1])
    parm <- c((sumObj$"estimates")[,1])    
    lenPV<-length(percVec)
#    strParm <- (unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2] 
#    strParm <- unique(obj[[9]][, ncol(obj[[9]]) - 1])  # second last column contains original curve levels
    strParm <- unique(obj$"data"[, ncol(obj$"data") - 1])  # second last column contains original curve levels
#    strParm <- strParm[apply(parmMat, 2, function(x){!any(is.na(x))})]
    

#    indexVec <- in1fct()
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
#
#    lenNV <- length(naVec[!is.na(naVec)])
#    if (lenNV>0)
#    {
#        parmMat <- parmMat[,-naVec[!is.na(naVec)]]
#        indexMat <- indexMat[,-naVec[!is.na(naVec)]]
#    }

    ncPM <- ncol(obj$"parmMat")  # [[10]])
    nrPM <- nrow(obj$"parmMat")  # [[10]])
    
    options(warn=-1)  # to avoid warnings when filling in matrix with elements in excess in the vector 
    indexMat <- t(matrix(NA, nrPM, ncPM))
    indexMat[!is.na(t(obj[[10]]))] <- 1:(nrPM*ncPM)
    indexMat <- t(indexMat)
    options(warn=0)
    
    naVec <- rep(FALSE, ncPM)
    for (i in 1:ncPM)
    {
        naVec[i] <- any(is.na(parmMat[, i]))
    }
    indexMat <- indexMat[, !naVec, drop=FALSE]
    parmMat <- parmMat[, !naVec, drop=FALSE] 
    strParm <- strParm[!naVec]

    ncPM2 <- ncol(parmMat)  # obj[[10]])
    nrPM2 <- nrow(parmMat)  # obj[[10]])
    indexMat <- matrix(1:(nrPM2*ncPM2), nrPM2, ncPM2, byrow = TRUE)   


    ## Finding out which parameter occurs most times; this determines the number of SI values
#    maxIndex <- 0
#    maxParm <- 0
#    for (i in indexVec)
#    {
#        PM <- parmMat[i,]
#        lenPM <- length(unique(PM))
#        if (lenPM > maxParm) {maxIndex <- match(unique(PM),PM); maxParm <- i}
#    }

#    nCol <- ncol(parmMat)
#    indexVec <- 1:nCol
#    for (i in 1:nCol)
#    {
#        if (any(is.na(parmMat[,i]))) {indexVec[i] <- NA}
#    }
#    indexVec <- indexVec[!is.na(indexVec)]
    indexVec <- 1:ncol(indexMat)    

#    lenEB <- maxParm
#    lenM <- length(indexMat[maxParm,])
    lenEB <- length(indexVec)
    lenM <- length(indexVec)


    ## Retrieving curve numbers 
#    parmName <- unique((unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2-1])[lenEB] 
#    compNamesTemp <- obj[[6]][grep(paste(parmName,":",sep=""), obj[[6]])]
#    compNames <- (unlist(strsplit(compNamesTemp, ":")))[(1:length(compNamesTemp))*2]
    compNames <- as.character(strParm)  # converting a factor


    ## Calculating SI values
#    numComp <- (lenPV*(lenPV+1)/2)*(lenM*(lenM-1)/2)
    numComp <- (lenPV*(lenPV-1)/2)*(lenM*(lenM-1)/2)
    siMat <- matrix(0, numComp, 4)
    matchVec <- rep(TRUE, numComp)
    dimNames <- rep("", numComp)
    
    degfree <- obj$"sumList"$"df.residual"  # df.residual(obj)  # obj$"summary"[6]  # sumObj[[4]][2]

    rowIndex <- 1
    for (i in 1:lenPV)
    {
        for (ii in 1:lenPV)
        {
            if (i>=ii) {next}
            pVec <- percVec[c(i,ii)]

            for (j in 1:lenM)
            {
                for (k in 1:lenM)
                {
                    if (j>=k) {next}
                    matchVec[rowIndex] <- (is.null(compMatch) || all(c(compNames[j],compNames[k])%in%compMatch))  # this is a canonical 2 as PAIRS are matched

                    jInd <- j
                    kInd <- k
                    if (reverse) {jInd <- k; kInd <- j}
                    
                    parmInd1 <- indexMat[, jInd]
                    parmInd2 <- indexMat[, kInd]
                    parmInd <- c(parmInd1, parmInd2)

                    varCov <- varMat[parmInd, parmInd]
                    parmChosen1 <- parmMat[, jInd]
                    parmChosen2 <- parmMat[, kInd]

#                    SIfctVal <- SIfct(c(in2fct(parmChosen1), in2fct(parmChosen2)), pVec)
                    SIeval <- SIlist(parmChosen1, parmChosen2, pVec, ...)
#                    print(SIeval)
                    siMat[rowIndex, 1] <- SIeval[[1]]
                    dimNames[rowIndex] <- paste(strParm[jInd], "/", strParm[kInd], ":", pVec[1], "/", pVec[2], sep="")

                    derEval1 <- SIeval[[2]]
                    derEval2 <- SIeval[[3]]
                    derEval <- c(derEval1, derEval2)

                    siMat[rowIndex, 2] <- sqrt(derEval%*%varCov%*%derEval)


                    ## Testing SI equal to 1
                    siMat[rowIndex, 3] <- (siMat[rowIndex, 1] - 1)/siMat[rowIndex, 2]
                    siMat[rowIndex, 4] <- (pt(-abs(siMat[rowIndex,3]), degfree))+(1-pt(abs(siMat[rowIndex, 3]), degfree))

                    rowIndex <- rowIndex+1
                }
            }
        }
    }
    dimnames(siMat)<-list(dimNames, c("Estimate", "Std. Error", "t-value", "p-value"))
    siMat <- siMat[matchVec, ,drop=FALSE]
 
    cat("\n") 
    printCoefmat(siMat)
    invisible(siMat)   
#    return(siMat[matchVec, ,drop=FALSE])
}
