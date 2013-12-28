"estfun.drc" <- function (x, ...) 
{
    if (identical(x$type, "binomial"))
    {  # not handling missing values
        nTotal <- weights(x)
        nObs <- x[["dataList"]][["resp"]] * nTotal
        fittedVal <- fitted(x)
        rval0 <- nObs/fittedVal + (nObs - nTotal)/(1 - fittedVal)
        rval <- x$deriv1 * rval0  
        rval
    }
    if (identical(x$type, "Poisson"))
    {  # not handling missing values
        resp <- x[["dataList"]][["resp"]]
        fittedVal <- fitted(x)
        rval0 <- resp / fittedVal - 1
        rval <- x$deriv1 * rval0  
        rval
    }    
    if (identical(x$type, "event"))
    {  # not handling missing values
        resp <- diff(c(x[["dataList"]][["resp"]], 1))
        fittedVal <- fitted(x)
        fittedVal[length(fittedVal)] <- 1
#        fittedVal2 <- c(fittedVal[-1], 1)  # assuming data ordered according to time
        rval0 <- resp / diff(fittedVal)
        lagDeriv1 <- apply(x$deriv1, 2, function(x){diff(x)})
        rval <- lagDeriv1 * rval0  
        rval
    }    
    if (identical(x$type, "continuous"))
    {    
        rval <- (weights(x) * residuals(x)) * x$deriv1
    }
    colnames(rval) <- names(coef(x))
    rval
}

"bread.drc" <- function (x, ...) 
{
#    if (identical(x$type, "binomial"))
#    {
#        breadMat <- vcov(x) * unlist(x$sumList[1])
#    }
    if (identical(x$type, "continuous"))
    {
#        breadMat <- summary(x)$cov.unscaled * unlist(x$sumList[1])
        breadMat <- vcov(x) / (summary(x)$rse[1]^2) * unlist(x$sumList[1])
    } else { ## Note: not checked for event time data!
        breadMat <- vcov(x) * unlist(x$sumList[1])
    }
    return(breadMat)
}

