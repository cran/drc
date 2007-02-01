"residuals.drc" <-
function(object, type = c("working", "standardised"), ...)
{
    type <- match.arg(type)
    
    if (type == "standardised")
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
