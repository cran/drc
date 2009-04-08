ED <- function (object, ...) UseMethod("ED", object)

"ED.drc" <-
function(object, respLev, interval = c("none", "delta", "fls", "tfls"), 
level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), 
type = c("relative", "absolute"), lref, uref, bound = TRUE, od = FALSE, display = TRUE, 
pool = TRUE, logBase = NULL, ...)
{
    interval <- match.arg(interval)
    reference <- match.arg(reference)
    type <- match.arg(type)
    
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
#    if ( (is.null(logBase)) && (!is.null(obj$"curve"[[2]])) )
#    {
#        logDose <- obj$"curve"[[2]]
#    }
#    if (inherits(object, "bindrc"))
#    {
#        ED2(object, percVec)
#    } else {

    ## Checking 'respLev' vector ... should be numbers between 0 and 100
    if ( (type == "relative") && (bound) ) 
    {
        if (any(respLev <= 0 | respLev >= 100)) 
        {
            stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
        }
    }
    lenPV <- length(respLev)

    ## Retrieving relevant quantities
    EDlist <- object$fct$"edfct"  
    if (is.null(EDlist)) {stop("ED values cannot be calculated")}      

#    assayNo <- object$data[, 3]  # obj[[9]][,3]
#    numAss <- length(unique(assayNo))

#    parmMat <- obj$"parmMat"  # [[10]]
#    sumObj <- summary(object, od = od)
#    varMat <- sumObj$"varMat"
#    vcMat <- vcov(object)
#    resVar <- sumObj$"resVar"  # [[1]]
#    varMat<-obj$"transformation"%*%sumObj$"varMat"%*%t(obj$"transformation")        
    
#    varMat <- obj[[12]]%*%sumObj[[2]]%*%t(obj[[12]])
#    varMat <- sumObj[[2]]
#    parm <- c((sumObj[[3]])[,1])
#    parm <- c((sumObj$"estimates")[,1])
#    strParm <- (unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2] 
#    strParm <- unique(obj[[9]][, ncol(obj[[9]]) - 1])  # second last column contains original curve levels

#    oData <- object$"data"
#    strParm <- unique(oData[, ncol(oData) - 1])  # second last column contains original curve levels
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
    
#    ncPM <- ncol(parmMat)  # obj[[10]])
#    nrPM <- nrow(parmMat)  # obj[[10]])
#    stop()
    
#    options(warn=-1)  # to avoid warnings when filling in matrix with elements in excess in the vector 
#    indexMat <- t(matrix(NA, nrPM, ncPM))
#    indexMat[!is.na(t(obj[[10]]))] <- 1:(nrPM*ncPM)
#    indexMat[!is.na(t(parmMat))] <- 1:(nrPM*ncPM)
#    indexMat <- t(indexMat)
#    options(warn=0)
#    print(indexMat)
#    stop()
    
#    naVec <- rep(FALSE, ncPM)
#    for (i in 1:ncPM)
#    {
#        naVec[i] <- any(is.na(parmMat[, i]))
#    }
##    indexMat <- indexMat[, !naVec, drop=FALSE]
#    parmMat <- parmMat[, !naVec, drop=FALSE] 
#    strParm <- strParm[!naVec]
##    print(indexMat)
#
#    ncPM2 <- ncol(parmMat)  # obj[[10]])
#    nrPM2 <- nrow(parmMat)  # obj[[10]])
#    indexMat <- matrix(1:(nrPM2*ncPM2), nrPM2, ncPM2, byrow = TRUE)   

#    indexMat0 <- object$"indexMat"
#    noNA <- complete.cases(t(indexMat0))
#    indexMat <- t((t(indexMat0))[noNA, , drop = FALSE])
#    parmMat0 <- object$"parmMat"  # [[10]]
#    parmMat <-  t((t(parmMat0))[noNA, , drop = FALSE])
      
#    strParm <- unique(object$dataList$curveid)
#    strParm <- strParm[noNA]
#    strParm <- colnames(parmMat0)    
    
    indexMat <- object$"indexMat"
    parmMat <- object$"parmMat"
    strParm <- colnames(parmMat)    
    vcMat <- vcov(object, od = od, pool = pool)

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


    ## Calculating ED values
    
    ## Defining vectors and matrices
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM    
    lenEB <- ncolIM    
    dimNames <- rep("", lenEB*lenPV)
    EDmat <- matrix(0, lenEB*lenPV, 2)
    oriMat <- matrix(0, lenEB*lenPV, 2)  
    
#    for (i in maxIndex)

    ## Skipping curve id if only one curve is present
    lenIV <- lenEB  # ncol(indexMat)
    if (length(unique(strParm)) == 1) 
    {
        strParm[1:lenIV] <- rep("", lenIV)
    } else {
        strParm <- paste(strParm, ":", sep = "")
    }

    rowIndex <- 1
    for (i in indexVec)
    {
        parmInd <- indexMat[, i]
        varCov <- vcMat[parmInd, parmInd]
        parmChosen <- parmMat[, i]

        for (j in 1:lenPV)
        {
            EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)
            
#            EDmat[rowIndex, 1] <- EDeval[[1]]
            EDval <- EDeval[[1]]
            dEDval <- EDeval[[2]]
            
            oriMat[rowIndex, 1] <- EDval
            oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
                   
            if (!is.null(logBase))
            {
                EDval <- logBase^(EDval)                
                dEDval <- EDval * log(logBase) * dEDval
            }
            EDmat[rowIndex, 1] <- EDval
            EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)

#            dimNames[rowIndex] <- paste(strParm[i], ":", percVec[j], sep="")
            dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
            rowIndex <- rowIndex + 1
        }
    }
    colNames <- c("Estimate", "Std. Error")
    
    ## Using t-distribution for continuous data
    ## (only under the normality assumption)
    if (object$"type" == "continuous")
    {
        qFct <- function(x) {qt(x, df.residual(object))}
    } else { # Otherwise the standard normal distribution is used
        qFct <- qnorm
    }

    if (interval == "delta")
    {
        ciMat <- matrix(0, lenEB*lenPV, 2)
        tquan <- qFct(1 - (1 - level)/2)        
#        ciMat[, 1] <- EDmat[, 1] - qnorm(level + (1-level)/2)*EDmat[, 2]
#        ciMat[, 2] <- EDmat[, 1] + qnorm(level + (1-level)/2)*EDmat[, 2]
        ciMat[, 1] <- EDmat[, 1] - tquan * EDmat[, 2]
        ciMat[, 2] <- EDmat[, 1] + tquan * EDmat[, 2]
        colNames <- c(colNames, "Lower", "Upper")
        ciLabel <- "Delta method"
    }
    if (interval == "tfls")
    {
        lsVal <- log(oriMat[, 1])
        lsdVal <- oriMat[, 2]/oriMat[, 1]
        tquan <- qFct(1 - (1 - level)/2)
                        
        ciMat <- matrix(0, lenEB*lenPV, 2)
        ciMat[, 1] <- exp(lsVal - tquan * lsdVal)
        ciMat[, 2] <- exp(lsVal + tquan * lsdVal)
        colNames <- c( colNames, "Lower", "Upper") 
        ciLabel <- "To and from log scale"       
    }
#    if ( (!is.null(logBase)) && (ci == "fls") )
    if (interval == "fls")
    {        
        ciMat <- matrix(0, lenEB*lenPV, 2)
        tquan <- qFct(1 - (1 - level)/2) 
        
        if (is.null(logBase)) 
        {
            logBase <- exp(1)
            EDmat[, 1] <- exp(EDmat[, 1])  # back-transforming log ED values
        }

#        oriVal <- log(EDeval[[1]], base = logBase)
#        oridVal1 <- EDeval[[2]]
#        oridVal2 <- sqrt(oridVal1%*%varCov%*%oridVal1)
        ciMat[, 1] <- logBase^(oriMat[, 1] - tquan * oriMat[, 2])
        ciMat[, 2] <- logBase^(oriMat[, 1] + tquan * oriMat[, 2])
        
        EDmat <- EDmat[, -2, drop = FALSE]  # standard errors not relevant        
        colNames <- c( colNames[-2], "Lower", "Upper")
        ciLabel <- "From log scale"  
    }
    if (!(interval == "none"))
    {
        EDmat <- as.matrix(cbind(EDmat, ciMat))
    } else {
        ciLabel <- NULL
    }   
    dimnames(EDmat) <- list(dimNames, colNames)
    resPrint(EDmat, "Estimated effective doses", interval, ciLabel, display = display)
}