"SI" <-
function(object, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, 
ci = c("none", "delta", "fieller", "fls"), level = ifelse(!(ci == "none"), 0.95, NULL), 
reference = c("upper", "control"), type = c("relative", "absolute"), logBase = NULL, ...)
{
     ## Matching the argument 'method'
     ci <- match.arg(ci)
     reference <- match.arg(reference)
     type <- match.arg(type)     
     
#     if( (ci == "fieller") && !(percVec == c(50, 50)) )
#     {
#         stop("Fieller's theorem only implemented for ED50_1/ED50_2")
#     }
     
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
#    SIlist <- obj$fct$"sifct"
#    if (is.null(SIlist)) {stop("SI values cannot be calculated")}

#    if (ci == "fls")
#    {
#        sifct <- createsifct(object$"fct"$"edfct")
#    } else {
#        sifct <- createsifct(object$"fct"$"edfct", logBase)
#    }

    if ( (is.null(logBase)) && (ci == "fls") )
    {
        stop("Argument 'logBase' not specified for ci = 'fls'")
    }
    sifct <- createsifct(object$"fct"$"edfct", logBase, ci == "fls")


    ## Checking contain of percVec vector ... should be numbers between 0 and 100
    if ( (type == "relative") && any(percVec<=0 | percVec>=100) ) 
    {
        stop("Percentages outside the interval [0, 100] not allowed")
    }

    if (missing(compMatch)) {matchNames <- FALSE} else {matchNames <- TRUE}


    ## Retrieving relevant quantities
    assayNo <- object$"data"[, 3]  # [[9]][,3]
    numAss <- length(unique(assayNo))

    sumObj<-summary(object, od = od)
    varMat <- sumObj$"varMat"
    
#    varMat<-obj$"transformation"%*%sumObj$"varMat"%*%t(obj$"transformation")    

#    varMat<-obj[[12]]%*%sumObj[[2]]%*%t(obj[[12]])
#    parm<-c((sumObj[[3]])[,1])
    parm <- c((sumObj$"estimates")[,1])    
    lenPV<-length(percVec)
#    strParm <- (unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2] 
#    strParm <- unique(obj[[9]][, ncol(obj[[9]]) - 1])  # second last column contains original curve levels
    strParm <- unique(object$"data"[, ncol(object$"data") - 1])  # second last column contains original curve levels
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

#    ncPM <- ncol(obj$"parmMat")  # [[10]])
#    nrPM <- nrow(obj$"parmMat")  # [[10]])
    
#    options(warn=-1)  # to avoid warnings when filling in matrix with elements in excess in the vector 
#    indexMat <- t(matrix(NA, nrPM, ncPM))
#    indexMat[!is.na(t(obj[[10]]))] <- 1:(nrPM*ncPM)
#    indexMat <- t(indexMat)
#    options(warn=0)
#    
#    naVec <- rep(FALSE, ncPM)
#    for (i in 1:ncPM)
#    {
#        naVec[i] <- any(is.na(parmMat[, i]))
#    }
#    indexMat <- indexMat[, !naVec, drop=FALSE]
#    parmMat <- parmMat[, !naVec, drop=FALSE] 
#    strParm <- strParm[!naVec]
#
#    ncPM2 <- ncol(parmMat)  # obj[[10]])
#    nrPM2 <- nrow(parmMat)  # obj[[10]])
#    indexMat <- matrix(1:(nrPM2*ncPM2), nrPM2, ncPM2, byrow = TRUE)   


    indexMat0 <- object$"indexMat"
    noNA <- complete.cases(t(indexMat0))
    indexMat <- t((t(indexMat0))[noNA, ])
    parmMat0 <- object$"parmMat"  # [[10]]
    parmMat <-  t((t(parmMat0))[noNA, ])
#    print(indexMat)


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
#    indexVec <- 1:ncol(indexMat)    

#    lenEB <- maxParm
#    lenM <- length(indexMat[maxParm,])

#    indexVec <- 1:ncol(indexMat)
    lenEB <- ncol(indexMat)  # length(indexVec)
    lenM <- ncol(indexMat)  # length(indexVec)


    ## Retrieving curve numbers 
#    parmName <- unique((unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2-1])[lenEB] 
#    compNamesTemp <- obj[[6]][grep(paste(parmName,":",sep=""), obj[[6]])]
#    compNames <- (unlist(strsplit(compNamesTemp, ":")))[(1:length(compNamesTemp))*2]
    compNames <- as.character(strParm)  # converting a factor


    ## Calculating SI values
#    numComp <- (lenPV*(lenPV+1)/2)*(lenM*(lenM-1)/2)
    numComp <- (lenPV*(lenPV-1)/2)*(lenM*(lenM-1)/2)
    
    if (!(ci == "none"))
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
    # obj$"sumList"$"df.residual"  # df.residual(obj)  # obj$"summary"[6]  # sumObj[[4]][2]


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
                    matchVec[rowIndex] <- (is.null(compMatch) || all(c(compNames[j],compNames[k])%in%compMatch))  
#                    # this is a canonical 2 as PAIRS are matched

                    jInd <- j
                    kInd <- k
                    if (reverse) {jInd <- k; kInd <- j}
                    
                    parmInd1 <- indexMat[, jInd]
                    parmInd2 <- indexMat[, kInd]
                    
                    splInd <- splitInd(parmInd1, parmInd2)

#                    parmInd <- c(parmInd1, parmInd2)
#                    varCov <- varMat[parmInd, parmInd]
                    
                    parmChosen1 <- parmMat[, jInd]
                    parmChosen2 <- parmMat[, kInd]
#                    parmChosen1 <- parmVec[parmInd1]  # parmVec is the vector parameter estimates
#                    parmChosen2 <- parmVec[parmInd2]

#                    print(parmChosen1) 
#                    print(parmChosen2)
#                    print(varMat)                    
#                    print(varCov)                    

#                    SIeval <- SIlist(parmChosen1, parmChosen2, pVec, ...)
                    SIeval <- 
                    sifct(parmChosen1, parmChosen2, pVec, 
                    splInd[[1]][, 1], splInd[[2]][, 1], splInd[[3]][, 1], splInd[[3]][, 2], reference, type, ...)
#                    print(SIeval)
                    indInOrder <- c(splInd[[1]][, 2], splInd[[2]][, 2], splInd[[3]][, 3])
                                       
                    SIval <- SIeval$"val"  # SIeval[[1]]
#                    derEval1 <- SIeval$"der1"  # SIeval[[2]]
#                    derEval2 <- SIeval$"der2"  # SIeval[[3]]
#                    derEval12 <- SIeval$"der12"  # SIeval[[4]]
#                    dSIval <- c(derEval1, derEval2, derEval12)                   
                    dSIval <- SIeval$"der"  # SIeval[[2]]
 
                    oriMat[rowIndex, 1] <- SIval
                    oriMat[rowIndex, 2] <- sqrt(dSIval%*%varMat[indInOrder, indInOrder]%*%dSIval)  # sqrt(dSIval%*%varCov%*%dSIval)

                    siMat[rowIndex, 1] <- SIval
                    rNames[rowIndex] <- paste(strParm[jInd], "/", strParm[kInd], ":", pVec[1], "/", pVec[2], sep="")

                    if (ci == "none")
                    {
#                        derEval1 <- SIeval[[2]]
#                        derEval2 <- SIeval[[3]]
#                        derEval <- c(derEval1, derEval2)
#                        siMat[rowIndex, 2] <- sqrt(derEval%*%varCov%*%derEval)
                        siMat[rowIndex, 2] <- oriMat[rowIndex, 2]  # sqrt(dSIval%*%varCov%*%dSIval)

                        ## Testing SI equal to 1
                        siMat[rowIndex, 3] <- (siMat[rowIndex, 1] - 1)/siMat[rowIndex, 2]
                        siMat[rowIndex, 4] <- (pt(-abs(siMat[rowIndex,3]), degfree))+(1-pt(abs(siMat[rowIndex, 3]), degfree))
                    }
                    if ( (ci == "delta") || (ci == "fls") )
                    {
#                        derEval1 <- SIeval[[2]]
#                        derEval2 <- SIeval[[3]]
#                        derEval <- c(derEval1, derEval2)
                        stErr <- oriMat[rowIndex, 2]  # sqrt(derEval%*%varCov%*%derEval)
                        tquan <- qt(1 - (1 - level)/2, degfree)  # df.residual(obj))
                        
                        siMat[rowIndex, 2] <- siMat[rowIndex, 1] - tquan * stErr
                        siMat[rowIndex, 3] <- siMat[rowIndex, 1] + tquan * stErr
                        ciLabel <- "Delta method"
                    }
                    if (ci == "tfls")
                    {
                        lsVal <- log(oriMat[rowIndex, 1])
                        lsdVal <- oriMat[rowIndex, 2]/oriMat[rowIndex, 1]
                        tquan <- qt(1 - (1 - level)/2, degfree)  # df.residual(obj))
                        
                        siMat[rowIndex, 2] <- exp(lsVal - tquan * lsdVal)
                        siMat[rowIndex, 3] <- exp(lsVal + tquan * lsdVal)
                        ciLabel <- "To and from log scale"
                    }
                    if ( (!is.null(logBase)) && (ci == "fls") )
                    {
                        siMat[rowIndex, 1] <- logBase^(siMat[rowIndex, 1])
                        siMat[rowIndex, 2] <- logBase^(siMat[rowIndex, 2])
                        siMat[rowIndex, 3] <- logBase^(siMat[rowIndex, 3])
                        ciLabel <- "From log scale"
                    }
                    if (ci == "fieller")
                    {
#                        lenp <- length(parmChosen1)  # model depedent!!!
#                        EDind <- c(parmInd1[lenp], parmInd2[lenp])
                        
                        vcMat <- matrix(NA, 2, 2)
                        vcMat[1, 1] <- SIeval$"der1"%*%varMat[parmInd1, parmInd1]%*%SIeval$"der1"
                        vcMat[2, 2] <- SIeval$"der2"%*%varMat[parmInd2, parmInd2]%*%SIeval$"der2"
                        vcMat[1, 2] <- SIeval$"der1"%*%varMat[parmInd1, parmInd2]%*%SIeval$"der2"
                        vcMat[2, 1] <- vcMat[1, 2]
#                        print(vcMat)
                        muVec <- c(SIeval$"valnum", SIeval$"valden")
                        
                        siMat[rowIndex, 2:3] <- fieller(muVec, degfree, vcMat, level = level)  
#                        fieller(c(parmChosen1[lenp], parmChosen2[lenp]), df.residual(obj), 
#                        varMat[EDind, EDind], level = level)                        
#                        varCov[c(lenp, 2*lenp), c(lenp, 2*lenp)], level = level)
                        ciLabel <- "Fieller"
                    }


                    rowIndex <- rowIndex+1
                }
            }
        }
    }
    dimnames(siMat)<-list(rNames, cNames)
    siMat <- siMat[matchVec, , drop = FALSE]
 
    cat("\n")
    cat("Estimated ratios of effect doses\n")
    if (!(ci == "none")) 
    {
        ciText <- paste("(", ciLabel, "-based confidence interval(s))\n", sep = "")
        cat(ciText)
    } 
    cat("\n")
    printCoefmat(siMat)
    invisible(siMat)   
#    return(siMat[matchVec, ,drop=FALSE])
}


"fieller" <-
function(mu, df, vcMat, level = 0.95, finney = FALSE, resVar)
{
    tper <- qt(1-0.5*(1-level), df)^2 
    
    if (!finney)
    {
        ## Based on the entry on Fieller's theorem 
        ##  in Encyclopedia of Statistical Sciences Vol. 3, p. 86 (1983)
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
