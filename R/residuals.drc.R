"residuals.drc" <-
function(object, ...)
{
    if (inherits(object, "bindrc"))
    {
        return(object$data$resp - fitted(object))
    } else {
    
        return(object[[7]][,2])
    }
}
