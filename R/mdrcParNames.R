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
        
    parmVec2 <- parmVec
    for (i in 1:length(parmVec))
    {
        pos <- regexpr("factor(collapse[, i])", parmVec[i], fixed=TRUE)
        if (pos>0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos+21), sep="")
        }
            
        pos <- regexpr("factor(assayNo)", parmVec[i], fixed=TRUE)
        if (pos>0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos+15), sep="")
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
