"residuals.drc" <-
function(object, typeRes = c("working", "standardised", "studentised"), ...)
{
    typeRes <- match.arg(typeRes)
    
    rawResiduals <- object$"dataList"$"resp" - fitted(object)
    
    if (identical(typeRes, "standardised"))
    {
        rstan <- object$"estMethod"$"rstanfct"
        if (is.null(rstan)) 
        {
            cat("Scale parameter fixed at 1. So working residuals are returned\n\n")
#            return(object$"predres"[, 2])
            return(rawResiduals)
            
#            stop("No standardisation available")
        } else {
#            return(object$"predres"[, 2] / rstan(object))  
            return(rawResiduals / rstan(object))  
            
#            return( object$"predres"[, 2] / sqrt(summary(object)$"resVar") )
        }
    } 
    
    if (identical(typeRes, "studentised"))
    {
        Xmat <- object$"deriv1"
        if (is.null(Xmat)) {stop("Studentised residuals not available")}
        
        Xprod <- solve(t(Xmat) %*% Xmat)
        nrowX <- nrow(Xmat)        
        Hdiag <- rep(NA, nrowX)
        for (i in 1:nrowX)
        {
            Hdiag[i] <- Xmat[i, ] %*% Xprod %*% t(Xmat[i, , drop = FALSE])
        }
#        print(length(Hdiag))
#        print(dim(object$"predres"))
        scaleEst0 <- summary(object)$"resVar"
        scaleEst <- ifelse(is.na(scaleEst0), 1, scaleEst0)  
        # to handle response types that are not continuous/quantitative
        
#        return(object$"predres"[, 2] / sqrt(scaleEst * (1 - Hdiag)))
        return(rawResiduals / sqrt(scaleEst * (1 - Hdiag)))
    }
    
    if (identical(typeRes, "working"))
    {
#        return(object$"predres"[, 2])
        return(rawResiduals)
    }
}

"scaleEst" <- function(object)
{



}
