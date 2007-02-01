"drmOpt" <- 
function(opfct, opdfct1, startVec, optMethod, constrained, warnVal, 
upperLimits, lowerLimits, errorMessage, maxIt, relTol, opdfct2 = NULL, parmVec) 
{
    ## Controlling the warnings
    options(warn = warnVal)
    
    
    ## Calculating hessian
    if (is.null(opdfct2)) {hes <- TRUE} else {hes <- FALSE}


    ## Optimising
    psVec <- abs(startVec)
    psVec[psVec < 1e-4] <- 1

    {if (!is.null(opdfct1))
    {
        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = "L-BFGS-B", 
            lower = lowerLimits, upper = upperLimits, control = list(maxit = maxIt)), silent = TRUE)
        } else {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = optMethod, 
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec)), silent = TRUE)
        }
        options(warn = 0)
        
        if (!inherits(nlsObj, "try-error")) 
        {
            nlsFit <- nlsObj
        } else {
#            stop("Convergence failed")
            warning("Convergence failed. The model was not fitted!", call. = FALSE)

            callDetail <- match.call()
            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}

            return(list(call = callDetail, parNames = parmVec, startVal = startVec))
        }
        if (!hes) {nlsFit$hessian <- opdfct2(parmVal)}

    } else {  # in case no derivatives are used

        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method="L-BFGS-B", 
            lower=lowerLimits, upper=upperLimits, 
            control=list(maxit=maxIt)), silent=TRUE)
        } else {
#            psVec <- abs(startVec)
#            psVec[psVec<1e-4] <- 1
            nlsObj <- try(optim(startVec, opfct, hessian = TRUE, method = optMethod, 
            control=list(maxit = maxIt, reltol = relTol, parscale = psVec)))  # , silent = TRUE)
#            nlsObj0 <- try(optim(startVec, opfct, method=optMethod, 
#            control=list(maxit=maxIt, reltol=relTol, parscale=psVec)), silent=TRUE)
#            nlsObj <- try(optim(nlsObj0$par, opfct, hessian=TRUE, method=optMethod, 
#            control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
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

    ## returning the fit
    return(nlsFit)
}
