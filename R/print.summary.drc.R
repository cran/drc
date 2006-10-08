"print.summary.drc" <-
function(x, ...)
{
    object <- x
    
    cat("\n")
    cat(paste("Model fitted: ", object$"fctName", "\n", sep=""))

    if (!is.null(object$"robust"))
    {
        cat("\n")
        cat("Robust estimation:", object$"robust", "\n")   
    }

    cat("\n")
    cat("Parameter estimates:\n\n")
    printCoefmat(object$"estimates")

    if ((!is.null(object$"resVar")) && (!(object$"type"=="binomial")))
    {
        cat("\n")
        cat("Estimated residual variance:", object$"resVar", "\n")  # , "\n\n")
    }

    if (!is.null(object$"varComp"))
    {
        cat("\n")
        cat("Estimated variance components:\n\n")
        printCoefmat(object$"varComp")
    }
    
    if (!is.null(object$"varParm"))  # summary of variance as power of mean
    {
        if (object$"varParm"$"type" == "varPower")
        {
            cat("\n")
            cat("Heterogeneity adjustment: variance is a power of mean\n\n")    
            if (dim(object$"varParm"$"estimates")[1] > 1)
            {
                printCoefmat(object$"varParm"$"estimates"[2, , drop=FALSE])  
                # only displaying power exponent, not residual variance
            } else {
                printCoefmat(object$"varParm"$"estimates")
            }
        }
        if (object$"varParm"$"type" == "hetvar")
        {
            cat("\n")
            cat("Estimated heterogeneous variances:\n\n")    
            printCoefmat(object$"varParm"$"estimates"[, 1, drop = FALSE])
        }        
    }
    
    if (!is.null(object$"boxcox"))  # summary of Box-Cox transformation
    {
        lambda <- object$"boxcox"[1]
        if (is.na(lambda)) 
        {
            # empty
        } else {
#            pVal <- format(object$"boxcox"[2], digits=3)
            boxcoxci <- c(format(object$"boxcox"[2], digits=3), format(object$"boxcox"[3], digits=3))

            cat("\n")
            cat("Heterogeneity adjustment: Box-Cox transformation\n\n")

            if (!is.na(boxcoxci[1]))
            {        
                cat("Estimated lambda:", lambda, "\n")
#                cat("P-value for test of null hypothesis that lambda=1:", pVal, "\n")
                ciStr <- paste("[", boxcoxci[1], ",", boxcoxci[2], "]", sep="")
                cat("Confidence interval for lambda:", ciStr, "\n\n")
            } else {
                cat("Specified lambda:", lambda, "\n\n")        
            }
        }
    }
         
    invisible(object)
}
