"idrm" <- function(x, y, curve, weights, fct, type)
{
    lenData <- length(y)
    dframe <- data.frame(y = y, x = x, w = weights)

#print(cbind(x,y,cid))

    ## Setting na.action option to default
    options(na.action = "na.omit")

    listData <- tapply(1:lenData, curve, function(t) {dframe[t, ]})

print(listData)

    laFct <- function(listDataElt, fct) {drm(y~x, weights=w, data=listDataElt, fct=fct, type=type)}  
    if (!is.list(fct[[1]]))
    {
#        fitList <- lapply(listData, multdrc, fct = fct, ...)
        fitList <- lapply(listData, laFct, fct)
    } else {
        fitList <- list()
        for (i in 1:length(fct))
        {
#            fitList[[i]] <- drm(listData[[i]], fct = fct[[i]], type = type, weights = weights)
            fitList[[i]] <- laFct(listData[[i]], fct[[i]])
        }
    }
    
#    ## The data set
#    if (!is.null(logDose)) 
#    {
#        dose <- origDose
#    }
#    dataSet <- data.frame(x, y, curve)  # , assayNoOld, weights)
#    names(dataSet) <- c(varNames, anName, anName, "weights")    
    
    retList <- list(fitList = fitList, fctList = fct, curveId = unique(curve))
    class(retList) <- c("idrm")
    return(retList)
#    return(fitList)
}

"summary.idrm" <- function(object)
{
    

}


"coef.idrm" <- function(object)
{
    lappFct <- function(t)
    {
        coefVec <- coef(t) 
        retVec <- c(coefVec, summary(t)$resVar)
        names(retVec) <- c(names(coefVec), "Res var")
        
        retVec
    }

#    coefList <- lapply(object$"fitList", function(t) {c(coef(t), summary(t)$resVar)})
    coefList <- lapply(object$"fitList", lappFct)

    if (!is.list(object$"fctList"[[1]]))
    {
        cl1 <- coefList[[1]]
        coefMat <- matrix(unlist(coefList), ncol = length(cl1), byrow = TRUE)
        colnames(coefMat) <- names(cl1)
        rownames(coefMat) <- object$"curveId"
    
        return(coefMat)
    } else {
        names(coefList) <- object$"curveId"
        return(coefList)
    }
}

