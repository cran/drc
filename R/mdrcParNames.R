"mdrcParNames" <- function(numNames, parNames, collapseList2)
{
    ## Retrieving names for parameters
    parmVecList <- list()
    for (i in 1:numNames)
    {
        colNames1 <- colnames(collapseList2[[i]])
        if (is.null(colNames1)) 
        {
            parmVecList[[i]] <- paste(parNames[i],"(Intercept)",sep=":")
        } else { 
            parmVecList[[i]] <- paste(parNames[i], colNames1, sep=":")
        }        
        parmVecList[[i]] <- (parmVecList[[i]])[1:ncol(collapseList2[[i]])]  # min(maxParm[i], length(colNames1))]
    }
    parmVec <- unlist(parmVecList)
#    print(parmVec)
        
    parmVec2 <- parmVec
    for (i in 1:length(parmVec))
    {
        posStr <- "factor(collapse[, i])"
        psLen <- nchar(posStr)
        pos <- regexpr(posStr, parmVec[i], fixed=TRUE)
        if (pos>0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos + psLen), sep = "")
        }
            
        posStr <- "factor(assayNo)"
        psLen <- nchar(posStr)
        pos <- regexpr(posStr, parmVec[i], fixed=TRUE)
        if (pos>0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos + psLen), sep = "")
        }

    }
#        parmVec <- parmVec2 

#    lenPV <- length(parmVec2)
#    parmVecA <- rep(0, lenPV)
#    parmVecB <- rep(0, lenPV)

#    splitList <- strsplit(parmVec2, ":")
#    for (i in 1:lenPV)
#    {
#        parmVecA[i] <- splitList[[i]][1]
#        parmVecB[i] <- splitList[[i]][2]    
#    }
    
    return(mdrcPNsplit(parmVec2, ":"))

#    return(list(parmVec2, parmVecA, parmVecB))
}

"mdrcPNsplit" <- function(parmVec, sep) 
{
    lenPV <- length(parmVec)
    parmVecA <- rep(0, lenPV)
    parmVecB <- rep(0, lenPV)

    splitList <- strsplit(parmVec, sep, fixed = TRUE)
    for (i in 1:lenPV)
    {
        parmVecA[i] <- splitList[[i]][1]
        
        lenSL <- length(splitList[[i]])
        parmVecB[i] <- paste(splitList[[i]][2:lenSL], collapse = "")  # 'paste' is needed in case several ":" occur
    }
    return(list(parmVec, parmVecA, parmVecB))
}
