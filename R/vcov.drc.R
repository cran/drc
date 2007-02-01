"vcov.drc" <-
function(object, ..., corr = FALSE)
{
    ## Retrieving the estimated variance-covariance matrix for the parameter estimates
    
    if (!corr)
    {
        summary(object)$"varMat"
    } else {
        vcMat <- summary(object)$"varMat"
        diage <- sqrt(diag(vcMat))
        scaMat <- outer(diage, diage)
        
        vcMat/scaMat
    }    
}
