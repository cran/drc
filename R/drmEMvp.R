"drmEMvp" <- function(dose, resp, multCurves)
{
    ## Defining the objective function
    lenData <- length(dose)

    opfct <- function(cVec)
    {
        lenc <- length(cVec)
        c1 <- cVec[-c(lenc-1,lenc)]
        c2 <- cVec[lenc-1]  # residual standard error
        c3 <- cVec[lenc]
                
        mc <- multCurves(dose, c1)
        return( 2*(lenData)*log(c2) + sum( c3*log(mc) + ((resp-mc)/(c2*(mc^(c3/2))))^2 ) )
    }

    
    ## Defining self starter function
    ssfct <- function(dataset)
    {
        xVal <- dataset[,1]
        yVal <- dataset[,2]

        yp <- log(tapply(yVal, xVal, var))
        xp <- log(tapply(yVal, xVal, mean))

        linModel <- lm(yp~xp)
        coefLM <- coef(linModel)
        retVal <- c(exp(coefLM[1]/2), coefLM[2])
#        if (retVal[1] < 0.001) {retVal[1] <- 0.001}  
#        # convergence fails if sigma becomes to small

        return(retVal)
    }


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        c( -(lenData/2)*log(2*pi)-object$"fit"$"value"/2, length(object$"fit"$"par") )
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance matrix and the fixed effects estimates
    rvfct <- function(object)
    {
        lp <- length(object$"fit"$"par")
        (object$"fit"$"par"[lp - 1])^2
    }

    vcovfct <- function(object)
    {
        solve(object$"fit"$"hessian")
    }
    
    parmfct <- function(fit, fixed=TRUE)
    {
        lp <- length(fit$"par")
        if (fixed) {fit$"par"[1:(lp - 2)]} else {fit$"par"}
    }

    rstanfct <- function(object)
    {
        sqrt(fitted(object)^(summary(object)$varParm$estimates[2,1])*summary(object)$"resVar")
    }

    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
                parmfct = parmfct, rstanfct = rstanfct))
}

"drmLOFvp" <- function()
{
    return(list(NULL, NULL))
}
