"colFct" <-
function(assayNo, noVec)
{
    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    }

    ## Collapsing a vector
    assayNo2 <- colConvert(assayNo)
    len <- length(assayNo2)
    assLev <- levels(factor(assayNo2))
    assayNo2[assayNo2%in%assLev[noVec]] <- paste(assLev[noVec],collapse="x")

    return(colConvert(assayNo2))
}
