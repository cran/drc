"summary.drc" <-
function(object, od = FALSE, ...)
{
#    if (inherits(object, "bindrc"))
#    {
#        sumbin(object, ...)
#    } else {


    ## Producing a summary of model fit
    sumVec1 <- object$fit  # object[[2]]
    sumVec2 <- object$summary  # object[[4]]
    parNames <- object$"parNames"[[1]]  # object[[6]]

    ## Calculating variance-covariance matrix from Hessian
    em <- object$"estMethod"
    parVec <- (em$parmfct)(object$fit, fixed = FALSE)
    notNA <- !is.na(parVec) 
    varMat <- (object$"scaleFct")( (em$vcovfct)(object) )
        
    ## Calculating estimated residual variance        
    if (!is.null(em$rvfct)) 
    {
        resVar <- (em$rvfct)(object)
#        resVar <- object$fit$value / df.residual(object)
    } else {
        resVar <- NULL
    }
    if (!is.null(resVar))
    {
        varMat.us <- varMat / (2*resVar)    
    } else {
        varMat.us <- NULL
    }
   
    if (od && (!is.null(object$"gofTest")) ) 
    { 
        varMat <- varMat*(object$"gofTest"[1]/object$"gofTest"[2])
    }
    estSE <- sqrt(diag(varMat))


    ## Calculating estimated standard errors for robust methods
    
    ## M-estimators
    if (!is.null(object$robust) && object$robust%in%c("metric trimming", "metric Winsorizing", "Tukey's biweight"))
    {
        psi.trimmed <- function(u, deriv = 0)
        {
            if (deriv == 0)
            {
                retVec <- u
                retVec[ abs(u) > 1.345 ] <- 0
            }
            if (deriv == 1)
            {
                retVec <- rep(1, length(u))
                retVec[ abs(u) > 1.345 ] <- 0
            }
            return(retVec)            
        }
    
    
        if (object$robust=="Tukey's biweight")
        {
            psifct <- psi.bisquare
        }
        if (object$robust=="metric Winsorizing")
        {
            psifct <- psi.huber
        }
        if (object$robust=="metric trimming")
        {
            psifct <- psi.trimmed
        }
    
        resVec <- residuals(object)
        psiprime <- psifct(resVec/sqrt(resVar), deriv=1)
        meanpp <- mean(psiprime)
        
        K <- 1 + length(parVec[notNA])*var(psiprime)/(object$summary[7]*meanpp^2)        
        w <- psifct(resVec/sqrt(resVar))
        s <- sum((resVec*w)^2)/object$summary[6]
        stddev <- sqrt(s)*(K/meanpp)       
        estSE <- sqrt(diag(solve(sumVec1$hessian[notNA, notNA]/mean(psiprime)*resVar)))*stddev
    }


    ## Forming a matrix of results        
    resultMat <- matrix(0, sum(notNA), 4, dimnames=list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))    
    resultMat[, 1] <- parVec[notNA]
    resultMat[, 2] <- estSE
    tempStat <- resultMat[, 1]/resultMat[, 2]
    resultMat[, 3] <- tempStat
    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (object$"type" == "continuous")
    {
        pFct <- function(x) {pt(x, df.residual(object))}
    } else {
        pFct <- pnorm
    }    
    resultMat[, 4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))

    ## Separating out variance parameters
    if (!is.null(object$"varParm"))
    {
        indexVec <- object$"varParm"$"index"
        varParm <- object$"varParm"

        estVec <- resultMat[-indexVec, , drop = FALSE]
        if (object$"varParm"$"type" == "varPower")
        {
            estVec[2, 3] <- (estVec[2, 1] - 0)/estVec[2, 2]  # testing the hypothesis theta=0
            estVec[2, 4] <- 2*pt(-abs(estVec[2, 3]), df.residual(object))
        }
        varParm$"estimates" <- estVec
        
        resultMat <- resultMat[indexVec,]
        varMat <- varMat[indexVec, indexVec]  # for use in ED/MAX/SI
    } else {
        varParm <- NULL
    }

    fctName <- deparse(object$call$fct)
    

    sumObj <- list(resVar, varMat, resultMat, object$"boxcox", fctName, object$"robust", varParm, object$"type", 
    df.residual(object), varMat.us, object$"fct"$"text", object$"fct"$"noParm")
    names(sumObj) <- c("resVar", "varMat", "coefficients", "boxcox", "fctName", "robust", "varParm", "type", 
    "df.residual", "cov.unscaled", "text", "noParm")
    class(sumObj) <- c("summary.drc")
    return(sumObj)
#    }
}
