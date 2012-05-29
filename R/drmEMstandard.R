"drmEMstandard" <- 
function(dose, resp, multCurves, weightVec, doseScaling = 1)
{

    ## Defining a helper function for calculating the variance-covariance matrix
#    vcFct <- function(beta0, beta, sigma2, len0)
#    {
#        vc <- (sigma2 / len0) * (beta %o% beta) / (beta0^4)
#        diag(vc) <- diag(vc) + sigma2 / (beta0^2)
#
#        return(vc)
#    }
    vcFct <- function(beta0, beta, len0)
    {
        vc <- (1 / len0) * (beta %o% beta) / (beta0^4)
        diag(vc) <- diag(vc) + (1 / (beta0^2))

        return(vc)
    }
    
    zeroDose <- dose < 1e-15  # hardcoded tolerance of 1e-15
    vcFct2 <- function(beta0, betaVec)
    {
        len0 <- weightVec[1]  # in case len0 is a vector
    
        vc <- (1 / len0) * (betaVec %o% betaVec) / (beta0^4)
        diag(vc) <- diag(vc) + (1 / (beta0^2))
        
#        zeroDose <- dose < doseTol
#        print(vc[!zeroDose, zeroDose])
#        print((1 / len0) * (-betaVec / (beta0^3)))
#        print(vc[!zeroDose, zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3)))
        
        vc[!zeroDose, zeroDose] <- vc[!zeroDose, zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3)) 
        vc[zeroDose, !zeroDose] <- vc[zeroDose, !zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3))
#        print(vc[zeroDose, zeroDose])
#        print(diag(vc[zeroDose, zeroDose]) + (1 / (len0 * beta0^2)) - (1 / (beta0^2)))
        diag(vc[zeroDose, zeroDose]) <- diag(vc[zeroDose, zeroDose]) + (1 / (len0 * beta0^2)) - (1 / (beta0^2))

        return(vc)
    }    


    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {
        f0 <- multCurves(0, c)[1]
        fVec <- multCurves(dose / doseScaling, c)
#        vcMat <- vcFct(f0, fVec, weightVec)   
        vcMat <- vcFct2(f0, fVec)   
        
        sum( (resp - fVec) %*% solve(vcMat) %*% (resp - fVec))
    }    

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )
        
        c(
        -object$"fit"$value + sum(log(gamma(resp+1))),
        object$"sumList"$"df.residual"
        )  # adding scale constant
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance matrix, and the parameter estimates
#    rvfct <- function(object)
#    {
#        object$"fit"$"value" / df.residual(object)  # object$"sumList"$"df.residual"
#    }
#
#    vcovfct <- function(object)
#    {
#        solve(object$fit$hessian)    
#    }
#

    # copied from drmEMls.R
    rvfct <- function(object)
    {
        object$"fit"$"value" / df.residual(object)
    }

    vcovfct <- function(object)
    {
        scaledH <- (object$"fit"$"hessian") / (2 * rvfct(object))
        invMat <- try(solve(scaledH), silent = TRUE)
    
        if (inherits(invMat, "try-error"))
        {
            ## More stable than 'solve' (suggested by Nicholas Lewin-Koh - 2007-02-12)
            ch <- try(chol(scaledH))
            if(inherits(ch, "try-error")) 
            {
                ch <- try(chol(0.99 * object$fit$hessian + 0.01 * diag(dim(object$fit$hessian)[1])))
            }
            ## Try regularizing if the varcov is unstable
            if(!inherits(ch, "try-error")) return(chol2inv(ch))
        } else {
            return(invMat)
        }
    } 
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFstandard" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
