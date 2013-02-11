"estfun.drc" <- function (x, ...) 
{
    rval <- residuals(x) * x$deriv1
    colnames(rval) <- names(coef(x))
    rval
}

"bread.drc" <- function (x, ...) 
{
    return(summary(x)$cov.unscaled * unlist(x$sumList[1]))
}

