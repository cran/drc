"mdrcPNsplit" <- function(parmVec, sep) 
{
    lenPV <- length(parmVec)
    parmVecA <- rep(0, lenPV)
    parmVecB <- rep(0, lenPV)

    splitList <- strsplit(parmVec, sep, fixed = TRUE)
    for (i in 1:lenPV)
    {
        parmVecA[i] <- splitList[[i]][1]
        parmVecB[i] <- splitList[[i]][2]
    }
    return(list(parmVec, parmVecA, parmVecB))
}
