"vcov.drc" <-
function(object, ...)
{
    ## Retrieving the estimated covariance matrix for the parameter estimates
    summary(object)$"varMat"
}
