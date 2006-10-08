"residuals.drc" <-
function(object, rstandard = FALSE, ...)
{
    if (rstandard)
    {
        rstan <- object$"estMethod"$"rstanfct"
        if (is.null(rstan)) 
        {
            stop("No standardisation available")
        } else {
            return(object$"predres"[, 2] / rstan(object))    
        }
    } else {
        return(object$"predres"[, 2])
    }
}
