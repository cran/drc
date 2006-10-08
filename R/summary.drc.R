"summary.drc" <-
function(object, od = FALSE, ...)
{
    if (inherits(object, "bindrc"))
    {
        sumbin(object, ...)
    } else {
    
#    if (inherits(object, "mixdrc"))
#    {
#        objectCopy <- object
#        class(objectCopy) <- c("nlme", "lme")
#        
#        sumObj <- summary(objectCopy)
#            
#        resultMat <- as.matrix(sumObj$tTable[, c(1,2,4,5)]) 
#
#        varComp <- matrix(as.numeric(VarCorr(objectCopy)[,1]))
#        colnames(varComp) <- "Variance"
#        rownames(varComp) <- rownames(VarCorr(objectCopy))
#
#        ll <- logLik(objectCopy)
#        loglik <- ll[1] 
#        degfre <- sumObj$dims$N - attr(ll, "df")
#
#        estimates <-  as.vector(sumObj$coefficients$fixed)
#        parNames <- rownames(resultMat)
#        varMat <- sumObj$varFix
#        
#        ## Defining return list
#        retList <- list(varComp, varMat, resultMat, c(loglik, degfre), parNames, "logistic with random effects")
#
#        names(retList) <- c("varComp", "varMat", "estimates", "loglik", "parNames", "class") 
#        class(retList) <- c("summary.drc")
#        return(retList)
#        
#    } else {
    

    ## Producing a summary of model fit

    sumVec1 <- object$fit  # object[[2]]
    sumVec2 <- object$summary  # object[[4]]
    parNames <- object$"parNames"[[1]]  # object[[6]]
    
#    parVec <- sumVec1$par


#    lenPV <- length(parmVec)
#    resultMat <- matrix(0, lenPV, 4, dimnames=list(parmVec, c("Estimate", "Std. Error", "t-value", "p-value")))
    
#    resultMat <- matrix(0, sum(notNA), 4, dimnames=list(parmVec, c("Estimate", "Std. Error", "t-value", "p-value")))    
#    resultMat[,1] <- sumVec1$par[notNA]

    ## Calculating estimated residual variance
#    resVar <- sumVec1$value/sumVec2[6]


    ## Calculating variance-covariance matrix from Hessian
#    notNA <- !is.na(sumVec1$par)    
#    notNA <- !is.na(parVec)    
#    varMat <- solve((sumVec1$hessian[notNA, notNA])*(1/resVar)/2)
#    estSE <- sqrt(diag(varMat))


    em <- object$"estMethod"
    parVec <- (em$parmfct)(object$fit, fixed=FALSE)
    notNA <- !is.na(parVec) 
        
    if (!is.null(em$rvfct)) 
    {
        resVar <- (em$rvfct)(object)
    } else {
        resVar <- NULL
    }
    varMat <- (object$"scaleFct")( (em$vcovfct)(object) )
   
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
#        print(K)
        
        w <- psifct(resVec/sqrt(resVar))
        s <- sum((resVec*w)^2)/object$summary[6]
        stddev <- sqrt(s)*(K/meanpp)
#        print(stddev)
        
        estSE <- sqrt(diag(solve(sumVec1$hessian[notNA, notNA]/mean(psiprime)*resVar)))*stddev
    }


    ## Forming a matrix of results        
    resultMat <- matrix(0, sum(notNA), 4, dimnames=list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))    
#    resultMat[,1] <- sumVec1$par[notNA]
    resultMat[,1] <- parVec[notNA]
    
    resultMat[,2] <- estSE

    resultMat[,3] <- resultMat[,1]/resultMat[,2]
    dfres <- object$"sumList"$"df.residual"
    resultMat[,4] <- (pt(-abs(resultMat[,3]), dfres))+(1-pt(abs(resultMat[,3]), dfres))


    ## Separating out variance parameters
    if (!is.null(object$"varParm"))
    {
        indexVec <- object$"varParm"$"index"
        varParm <- object$"varParm"
#        varParm$"estimates" <- resultMat[-indexVec, , drop=FALSE]

        estVec <- resultMat[-indexVec, , drop=FALSE]
        if (object$"varParm"$"type" == "varPower")
        {
            estVec[2, 3] <- (estVec[2, 1] - 0)/estVec[2, 2]  # testing the hypothesis theta=0
            estVec[2, 4] <- 2*pt(-abs(estVec[2, 3]), dfres)
        }
        varParm$"estimates" <- estVec
        
        resultMat <- resultMat[indexVec,]
        varMat <- varMat[indexVec, indexVec]  # for use in ED/MAX/SI
    } else {
        varParm <- NULL
    }


#    ## Calculating the value of the log-likelihood
#    if (object$"type"=="continuous")
#    {
##        degfre <- sumVec2[6]
##        loglik <- (-(degfre/2)*(log(2*pi)+log(sumVec2[5])-log(degfre)+1))
#        
#        loglik <- (em$llfct)(object)
#    }
#    if (object$"type"=="binomial")
#    {
#        degfre <- sumVec2[6]
#        
#        total <- (object$"data")[,5]
#        success <- total*(object$"data")[,2]    
#        loglik <- sum(log(choose(total, success))) - object$fit$ofvalue
#    }


    fctName <- deparse(object$call$fct)

#    classObj <- class(object)
#    if (inherits(object, "list")) 
#    {
#        modelClass <- classObj[length(classObj)-1]
#    } else {
#        modelClass <- classObj[length(classObj)]
#    }

    sumObj <- list(resVar, varMat, resultMat, object$"boxcox", fctName, object$"robust", varParm, object$"type")
    names(sumObj) <- c("resVar", "varMat", "estimates", "boxcox", "fctName", "robust", "varParm", "type")
    class(sumObj) <- c("summary.drc")
    return(sumObj)
    }
#    }
}
