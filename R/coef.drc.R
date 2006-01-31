"coef.drc" <-
function(object, ...)
{
    ## Retrieving the parameter estimates
    retVec <- object$fit$par
    names(retVec) <- object$parNames[[1]]

    return(retVec)
}
