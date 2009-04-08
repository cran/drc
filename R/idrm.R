"idrm" <- function(x, y, curveid, weights, fct, type)
{
    oneFunction <- !is.list(fct[[1]])
 
    ## Fitting models for each curve
    fitList <- list() 
    uniCur <- unique(curveid)   
    numCur <- length(uniCur)
    for (i in 1:numCur)
    {
        if (oneFunction)
        {
            fitList[[i]] <- drm(y~x, curveid, weights = weights, fct = fct, type = type, 
            subset = curveid == uniCur[i], separate = FALSE)
        } else {
            tempFitlist <- list()
            for (j in 1:length(fct))
            {
                tempFitlist[[j]] <- drm(y~x, curveid, weights = weights, fct = fct[[j]], 
                type = type, subset = curveid == uniCur[i], separate = FALSE)
            }
            fitList[[i]] <- tempFitlist
        }
    }

    retList <- fitList[[1]]
    if (oneFunction)
    {
        parmMat <- fitList[[1]]$"parmMat"
        nlsFit <- fitList[[1]]$"nlsFit"
        
        parNames <- fitList[[1]]$"parNames"
        numPar <- length(parNames[[1]])        
        parNames[[3]] <- rep(uniCur[i], numPar)
               
        if (numCur > 1)
        {
            nlsFit[[1]] <- nlsFit
            pnList <- list()
            pnList[[1]] <- parNames
            pnList[[1]][[3]] <- rep(uniCur[1], numPar)
            
            
            for (i in 2:numCur)
            {
                parmMat <- cbind(parmMat, fitList[[i]]$"parmMat")
                nlsFit[[i]] <- fitList[[i]]$"nlsFit"
                pnList[[i]] <- fitList[[i]]$"parNames"
                pnList[[i]][[3]] <- rep(uniCur[i], numPar)
            }
        }
        retList$"parmMat" <- parmMat 

        bVec <- as.vector(unlist(lapply(pnList, function(x){x[[2]]})))
        cVec <- as.vector(unlist(lapply(pnList, function(x){x[[3]]}))) 
        aVec <- paste(bVec, cVec, sep = ":")
        retList$"parNames" <- list(aVec, bVec, cVec)       
    
        retList$"indexMat" <- matrix(c(1:(numCur * numPar)), numPar, numCur)
        names(fitList) <- uniCur
        retList$"objList" <- fitList
    
        coefVec <- as.vector(unlist(lapply(fitList, function(x){x$"fit"$"par"})))
        names(coefVec) <- aVec
        retList$"coefficients" <- coefVec
        
        retList$"df.residual" <- sum(unlist(lapply(fitList, function(x){x$"df.residual"})))
        retList$"minval" <- sum(unlist(lapply(fitList, function(x){x$"fit"$"value"})))
        
        retList$"fit" <- nlsFit
    
        
    } else {
    
    }
        
    class(retList) <- c("drc")
    return(retList)
}

#
#"summary.idrm" <- function(object)
#{
#    
#
#}
#
#
#"coef.idrm" <- function(object)
#{
#    lappFct <- function(t)
#    {
#        coefVec <- coef(t) 
#        retVec <- c(coefVec, summary(t)$resVar)
#        names(retVec) <- c(names(coefVec), "Res var")
#        
#        retVec
#    }
#
##    coefList <- lapply(object$"fitList", function(t) {c(coef(t), summary(t)$resVar)})
#    coefList <- lapply(object$"fitList", lappFct)
#
#    if (!is.list(object$"fctList"[[1]]))
#    {
#        cl1 <- coefList[[1]]
#        coefMat <- matrix(unlist(coefList), ncol = length(cl1), byrow = TRUE)
#        colnames(coefMat) <- names(cl1)
#        rownames(coefMat) <- object$"curveId"
#    
#        return(coefMat)
#    } else {
#        names(coefList) <- object$"curveId"
#        return(coefList)
#    }
#}

