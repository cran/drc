"PR" <- function(object, xVec, ...)
{
    lenXV <- length(xVec)
    
    curveId <- as.character(unique(object$data[, 3]))
    lenCI <- length(curveId)

    retMat <- predict(object, data.frame(xVec, rep(curveId, rep(lenXV, lenCI))), ...)
    
    if (lenCI > 1)
    {
        rownames(retMat) <- paste(rep(curveId, rep(lenXV, lenCI)), rep(as.character(xVec), lenCI), sep = ":")
    } else {
        rownames(retMat) <- rep(as.character(xVec), lenCI)
    }
    
    return(retMat)
}
