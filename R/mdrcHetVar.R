"mdrcHetVar" <- function(dose, resp, multCurves, vvar, anovaYes = FALSE)
{
    ## Defining goodness-of-fit function
    anovaTest <- contAnovaTest()
    gofTest <- NULL
    if (anovaYes) {return(list(anovaTest = anovaTest, gofTest = gofTest))}


    ## Defining the objective function
    lenData <- length(dose)
    lenvvar <- length(unique(vvar))
    repVar <- table(vvar)[unique(as.character(vvar))]

    opfct <- function(cVec)
    {
        lenc <- length(cVec)
        c1 <- cVec[1:(lenc-lenvvar)]
        c2 <- cVec[(lenc-lenvvar+1):lenc] 
        varVec <- rep(c2, repVar)        
                
        mc <- multCurves(dose, c1)
        return( sum( log(varVec) + ((resp-mc)^2)/varVec ))
    }

    
    ## Defining self starter function
    ssfct <- function(dataset)
    {
#        xVal <- dataset[,1]
        yVal <- dataset[,2]

        return( tapply(yVal, vvar, var)/lenData )
    }


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        c( -(lenData/2)*log(2*pi)-object$"fit"$"value"/2, length(object$"fit"$"par") )
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance matrix and the fixed effects estimates
    rvfct <- function(object)
    {
        NULL
    }

    vcovfct <- function(object)
    {
        solve(object$"fit"$"hessian")
    }
    
    parmfct <- function(fit, fixed=TRUE)
    {
        lp <- length(fit$"par")
        if (fixed) {fit$"par"[1:(lp - lenvvar)]} else {fit$"par"}
    }


    return(list(llfct=llfct, opfct=opfct, ssfct=ssfct, rvfct=rvfct, vcovfct=vcovfct, parmfct=parmfct))
}
