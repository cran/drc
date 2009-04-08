"SI" <-
function(object, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, 
interval = c("none", "delta", "fieller", "fls"), level = ifelse(!(interval == "none"), 0.95, NULL), 
reference = c("control", "upper"), type = c("relative", "absolute"),
display = TRUE, pool = TRUE, logBase = NULL, ...)
{
     ## Matching the argument 'method'
     interval <- match.arg(interval)
     reference <- match.arg(reference)
     type <- match.arg(type)     

    if ( (is.null(logBase)) && (interval == "fls") )
    {
        stop("Argument 'logBase' not specified for interval = 'fls'")
    }
    sifct <- createsifct(object$"fct"$"edfct", logBase, identical(interval, "fls"))

    ## Checking contain of percVec vector ... should be numbers between 0 and 100
    if ( (type == "relative") && any(percVec<=0 | percVec>=100) ) 
    {
        stop("Percentages outside the interval [0, 100] not allowed")
    }

    if (missing(compMatch)) {matchNames <- FALSE} else {matchNames <- TRUE}

    lenPV <- length(percVec)

    ## Retrieving relevant quantities
    indexMat <- object$"indexMat"
    lenEB <- ncol(indexMat)    
    parmMat <- object$"parmMat"
    strParm <- colnames(parmMat)    
    varMat <- vcov(object, od = od, pool = pool)
#    compNames <- as.character(strParm)  # converting a factor

    ## Calculating SI values
    numComp <- (lenPV*(lenPV-1)/2)*(lenEB * (lenEB - 1) / 2)
    
    if (!identical(interval, "none"))
    {
        siMat <- matrix(0, numComp, 3)
        cNames <- c("Estimate", "Lower", "Upper")
    
    } else {
        siMat <- matrix(0, numComp, 4)
        cNames <- c("Estimate", "Std. Error", "t-value", "p-value")
    }
    matchVec <- rep(TRUE, numComp)
    rNames <- rep("", numComp)
    oriMat <- matrix(0, numComp, 2)    
    degfree <- df.residual(object)  
    rowIndex <- 1
    for (i in 1:lenPV)
    {
        for (ii in 1:lenPV)
        {
            if (i>=ii) {next}
            pVec <- percVec[c(i, ii)]

            for (j in 1:lenEB)
            {
                for (k in 1:lenEB)
                {
                    if (j>=k) {next}
                    matchVec[rowIndex] <- (is.null(compMatch) || all(c(strParm[j], strParm[k]) %in% compMatch))  
 
                    jInd <- j
                    kInd <- k
                    if (reverse) 
                    {
                        jInd <- k; kInd <- j; pVec <- pVec[c(2, 1)]
                    }
                   
                    parmInd1 <- indexMat[, jInd]
                    parmInd2 <- indexMat[, kInd]
                    
                    splInd <- splitInd(parmInd1, parmInd2)
                    
                    parmChosen1 <- parmMat[, jInd]
                    parmChosen2 <- parmMat[, kInd]

                    SIeval <- 
                    sifct(parmChosen1, parmChosen2, pVec, 
                    splInd[[1]][, 1], splInd[[2]][, 1], splInd[[3]][, 1], splInd[[3]][, 2], reference, type, ...)

                    indInOrder <- c(splInd[[1]][, 2], splInd[[2]][, 2], splInd[[3]][, 3])
                                       
                    SIval <- SIeval$"val"  # SIeval[[1]]
                    dSIval <- SIeval$"der"  # SIeval[[2]]
 
                    oriMat[rowIndex, 1] <- SIval
                    oriMat[rowIndex, 2] <- sqrt(dSIval%*%varMat[indInOrder, indInOrder]%*%dSIval)  # sqrt(dSIval%*%varCov%*%dSIval)

                    siMat[rowIndex, 1] <- SIval
                    rNames[rowIndex] <- paste(strParm[jInd], "/", strParm[kInd], ":", pVec[1], "/", pVec[2], sep="")

                    ## Using t-distribution for continuous data
                    ##  only under the normality assumption
                    if (identical(object$"type", "continuous"))
                    {
                        qFct <- function(x) {qt(x, degfree)}
                        pFct <- function(x) {pt(x, degfree)}
                    } else {
                        qFct <- qnorm
                        pFct <- pnorm
                    }

                    if (identical(interval, "none"))
                    {
                        siMat[rowIndex, 2] <- oriMat[rowIndex, 2]  # sqrt(dSIval%*%varCov%*%dSIval)

                        ## Testing SI equal to 1
                        tempStat <- (siMat[rowIndex, 1] - 1)/siMat[rowIndex, 2]
                        siMat[rowIndex, 3] <- tempStat
                        siMat[rowIndex, 4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))
                    }
                    if ( (identical(interval, "delta")) || (identical(interval, "fls")) )
                    {
                        stErr <- oriMat[rowIndex, 2]  # sqrt(derEval%*%varCov%*%derEval)
                        tquan <- qFct(1 - (1 - level)/2)
                        
                        siMat[rowIndex, 2] <- siMat[rowIndex, 1] - tquan * stErr
                        siMat[rowIndex, 3] <- siMat[rowIndex, 1] + tquan * stErr
                        ciLabel <- "Delta method"
                    }
                    if (identical(interval, "tfls"))
                    {
                        lsVal <- log(oriMat[rowIndex, 1])
                        lsdVal <- oriMat[rowIndex, 2]/oriMat[rowIndex, 1]
                        tquan <- qFct(1 - (1 - level)/2)
                        
                        siMat[rowIndex, 2] <- exp(lsVal - tquan * lsdVal)
                        siMat[rowIndex, 3] <- exp(lsVal + tquan * lsdVal)
                        ciLabel <- "To and from log scale"
                    }
                    if ((!is.null(logBase)) && (identical(interval, "fls")))
                    {
                        siMat[rowIndex, 1] <- logBase^(siMat[rowIndex, 1])
                        siMat[rowIndex, 2] <- logBase^(siMat[rowIndex, 2])
                        siMat[rowIndex, 3] <- logBase^(siMat[rowIndex, 3])
                        ciLabel <- "From log scale"
                    }
                    if (identical(interval, "fieller"))  # using t-distribution
                    {
                        vcMat <- matrix(NA, 2, 2)
                        vcMat[1, 1] <- SIeval$"der1"%*%varMat[parmInd1, parmInd1]%*%SIeval$"der1"
                        vcMat[2, 2] <- SIeval$"der2"%*%varMat[parmInd2, parmInd2]%*%SIeval$"der2"
                        vcMat[1, 2] <- SIeval$"der1"%*%varMat[parmInd1, parmInd2]%*%SIeval$"der2"
                        vcMat[2, 1] <- vcMat[1, 2]
                        muVec <- c(SIeval$"valnum", SIeval$"valden")
                        
                        siMat[rowIndex, 2:3] <- fieller(muVec, degfree, vcMat, level = level)  
                        ciLabel <- "Fieller"
                    }


                    rowIndex <- rowIndex+1
                }
            }
        }
    }
    dimnames(siMat) <- list(rNames, cNames)
    siMat <- siMat[matchVec, , drop = FALSE]

    resPrint(siMat, "Estimated ratios of effect doses\n", interval, ciLabel, display = display)
    
#    if (display)
#    { 
#        cat("\n")
#        cat("Estimated ratios of effect doses\n")
#        if (!(ci == "none")) 
#        {
#            ciText <- paste("(", ciLabel, "-based confidence interval(s))\n", sep = "")
#            cat(ciText)
#        } 
#        cat("\n")
#        printCoefmat(siMat)
#    }
#    invisible(siMat)   
}


"fieller" <-
function(mu, df, vcMat, level = 0.95, finney = FALSE, resVar)
{
    tper <- qt(1-0.5*(1-level), df)^2 
    
    if (!finney)
    {
        ## Based on the entry on Fieller's theorem 
        ##  in Encyclopedia of Statistical Sciences Vol. 3 (1983), p. 86 
        ##  essentially same formula as in Finney (see below)

        mup <- prod(mu)
    
        fVec0 <- mup - tper*vcMat[1,2]
        y2 <- mu[2]^2    
        fVec <- (fVec0)^2 - (mu[1]^2 - tper*vcMat[1,1])*(y2 - tper*vcMat[2,2])
    
        denom <- y2 - tper*vcMat[2,2]
        lowerL <- (fVec0 - sqrt(fVec))/denom
        upperL <- (fVec0 + sqrt(fVec))/denom
        
    } else {
    
        ## Using the formula 
        ##  in Finney: Statistical Method in Biological Assay p. 81 (3rd edition, 1978)
        ##  OOPS: uses the estimated residual variance
        fac <- sqrt(tper)*sqrt(resVar)/mu[2]
        g <- tper*vcMat[2,2]/(mu[2]^2)
        if (g >= 1) {stop("Fieller's theorem not useful!")} 
        ratio <- mu[1]/mu[2]
    
        v11 <- vcMat[1,1]/(resVar)
        v12 <- vcMat[1,2]/(resVar)
        v22 <- vcMat[2,2]/(resVar)
        innerBr <- g*(v11 - (v12^2)/v22)
        inBr <- v11 - 2*ratio*v12 + (ratio^2)*v22 - innerBr

        firstTerm <- ratio - g*vcMat[1,2]/vcMat[2,2]
        secondTerm <- fac*sqrt(inBr)
        denom <- 1 - g
        lowerL <- (firstTerm - secondTerm)/denom
        upperL <- (firstTerm + secondTerm)/denom
    }
    return(c(lowerL, upperL))
}

"splitInd"  <- function(ind1, ind2)
{
    matchVec1 <- ind1%in%ind2
    matchVec2 <- ind2%in%ind1
    
#    inCommon <- list(pos1 = (1:length(ind1))[matchVec1], pos2 = (1:length(ind2))[matchVec2], val = ind1[matchVec1])
#
    lmv1 <- sum(matchVec1)
    if (lmv1 > 0.01)
    {
        inCommon <- matrix( c( (1:length(ind1))[matchVec1], (1:length(ind2))[matchVec2], ind1[matchVec1]), lmv1, 3)
    } else {
        inCommon <- NULL
    }
    
#    only1 <- list(pos = (1:length(ind1))[!matchVec1], val = ind1[!matchVec1])
    only1 <- matrix( c( (1:length(ind1))[!matchVec1], ind1[!matchVec1] ), sum(!matchVec1), 2)
    
#    only2 <- list(pos = (1:length(ind2))[!matchVec2], val = ind2[!matchVec2])
    only2 <- matrix( c( (1:length(ind2))[!matchVec2], ind2[!matchVec2] ), sum(!matchVec2), 2)

    return(list(only1, only2, inCommon))
}

createsifct <- function(edfct, logBase = NULL, fls = FALSE)
{
    if (is.null(edfct)) 
    {
        stop("SI values cannot be calculated")
    } else {
        
        if (!fls)
        {
            if (is.null(logBase))
            {
                "sifct" <- function(parm1, parm2, pair, ind1, ind2, cmonInd1, cmonInd2, reference, type, ...)
                {
                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- ED1[[2]]
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- ED2[[2]]
        
                    SIpair <- ED1v/ED2v  # calculating the SI value
        
                    SIder1 <- ED1d/ED2v
                    SIder2 <- (-ED2d/ED2v)*SIpair
                    SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)
#                    SIder12 <- ED1d/ED2v - (ED2d/ED2v)*SIpair

                    return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),
                    der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
                }
            } else {
        
                "sifct" <- function(parm1, parm2, pair, ind1, ind2, cmonInd1, cmonInd2, reference, type, ...)
                {
                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- ED1[[2]]
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- ED2[[2]]
        
                    SIpair <- logBase^(ED1v - ED2v)  # calculating the SI value
        
                    SIder1 <- SIpair*log(logBase)*ED1d
                    SIder2 <- SIpair*log(logBase)*(-ED2d)
                    SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)

                    return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),
                    der1 = (log(logBase)*logBase^ED1v)*ED1d, der2 = (log(logBase)*logBase^ED2v)*ED2d, 
                    valnum = logBase^ED1v, valden = logBase^ED2v))
                }        
            }
        } else {
            
            "sifct" <- function(parm1, parm2, pair, ind1, ind2, cmonInd1, cmonInd2, reference, type, ...)
            {
                ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                ED1v <- ED1[[1]]
                ED1d <- ED1[[2]]
        
                ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                ED2v <- ED2[[1]]
                ED2d <- ED2[[2]]
        
            
                SIpair <- ED1v - ED2v  # calculating the log SI value
        
                SIder1 <- ED1d
                SIder2 <- -ED2d
                SIder12 <- commonParm(SIder1, SIder2, cmonInd1, cmonInd2)

                return(list(val = SIpair, der = c(SIder1[ind1], SIder2[ind2], SIder12),
                der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
            }
        }        
        return(sifct)
    }
}

commonParm <- function(SIder1, SIder2, cmonInd1, cmonInd2)
{
    lind1 <- length(cmonInd1)
    retVec <- rep(NA, lind1)
    for (i in 1:lind1)
    {
        retVec[i] <- SIder1[cmonInd1[i]] + SIder2[cmonInd2[i]]
    
    }
    return(retVec)
}

commatFct <- function(object, compMatch)
{
    parmMat <- object$parmMat

    if (!is.null(compMatch))
    {
        return(parmMat[, (colnames(parmMat) %in% c(compMatch[1], compMatch[2])), drop = FALSE ])
    } else {
        parmMat
    }
}
