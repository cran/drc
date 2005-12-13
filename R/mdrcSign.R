"mdrcSign" <- function(dose, resp, assayNo = rep(1, length(dose)) )
{

    lenUA <- length(unique(assayNo))
    diffVec <- rep(0, lenUA)
    for (i in 1:lenUA)
    {
        indVec <- assayNo==i

        minDose <- (dose[indVec])[which.min(resp[indVec])]
        maxDose <- (dose[indVec])[which.max(resp[indVec])]
        
        diffVec[i] <- maxDose - minDose
    }
    
    ## Majority vote: returns 1 if curves mostly are increasing and -1 if most decreasing
    if (sum(diffVec>0) > (lenUA/2 - (1e-12)) ) {return(1)} else {return(-1)}
}