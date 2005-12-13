"mdrcOpt" <- function(opfct, startVec, optMethod, derFlag, constrained, warnVal, upperLimits, lowerLimits, errorMessage, maxIt, relTol) 
{

    ## Optimising
    options(warn = warnVal)
    {if (derFlag)
    {
        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, opfctDer, method="L-BFGS-B", lower=lowerLimits, upper=upperLimits, control=list(maxit=maxIt)), silent=TRUE)
        } else {
            nlsObj <- try(optim(startVec, opfct, opfctDer, method=optMethod, control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
        }
        options(warn = 0)
        
        if (!inherits(nlsObj, "try-error")) 
        {
            nlsFit <- nlsObj
        } else {
#            stop("Convergence failed")
            warning("Convergence failed. The model was not fitted!", call.=FALSE)

            callDetail <- match.call()
            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}

            return(list(call = callDetail, parNames = parmVec, startVal = startVec))
        }
        nlsFit$hessian <- opfctDer2(parmVal)

    } else {  # in case no derivatives are used

        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method="L-BFGS-B", lower=lowerLimits, upper=upperLimits, 
                                control=list(maxit=maxIt)), silent=TRUE)
        } else {
            psVec <- abs(startVec)
            psVec[psVec<1e-4] <- 1
            nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method=optMethod, control=list(maxit=maxIt, reltol=relTol, parscale=psVec)), silent=TRUE)
#            nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method=optMethod, control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
        }
        options(warn = 0)
        
        if (!inherits(nlsObj, "try-error")) 
        {
            nlsFit <- nlsObj
        } else {  # to avoid an error if used in a loop
            if (errorMessage) {stop("Convergence failed")} else {warning("Convergence failed. The model was not fitted!", call. = FALSE)}

            callDetail <- match.call()
            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}
            return(list(call=callDetail, parNames=parmVec, startVal=startVec))
        }
    }}
    nlsFit$ofvalue <- nlsFit$value
#    options(warn=0)

    return(nlsFit)
}
