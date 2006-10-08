"mdrcLs" <- function(dose, resp, multCurves, startVec, robustFct, weights, rmNA, anovaYes = FALSE)
{
    ## Defining goodness-of-fit function
    anovaTest <- contAnovaTest()
    gofTest <- NULL
    if (anovaYes) {return(list(anovaTest = anovaTest, gofTest = gofTest))}


    ## Defining the objective function
    opfct <- function(c)
    {
        sum( robustFct( (resp - multCurves(dose, c))*weights ), na.rm = rmNA)  # weights enter multiplicatively!
    }

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        degfre <- object$"sumList"$"lenData"  # "df.residual"  # object$summary[6]
        c( -(degfre/2)*(log(2*pi)+log(object$"fit"$"value")-log(degfre)+1), length(object$"fit"$"par") + 1 )
    }
    
    
    ## Defining functions returning the residual variance, the variance-covariance matrix and the fixed effects estimates
    rvfct <- function(object)
    {
        object$"fit"$"value"/object$"sumList"$"df.residual"
    }

    vcovfct <- function(object)
    {
        solve((object$"fit"$"hessian")*(1/rvfct(object))/2)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }

    rstanfct <- function(object)
    {
        rep(1, object$"sumList"$"lenData")*sqrt(summary(object)$"resVar")
    }


    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct,
                parmfct = parmfct, rstanfct = rstanfct))
}
