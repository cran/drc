"drmEMbinomial" <- 
function(dose, resp, multCurves, startVec, robustFct, weights, rmNA, zeroTol = 1e-12, doseScaling = 1)
{
    ## Finding indices for doses that give contribution to likelihood function
    iv <- ( (multCurves(dose/doseScaling, startVec) > zeroTol) & (multCurves(dose/doseScaling, startVec) < 1-zeroTol) )


    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
        prob <- (multCurves(dose / doseScaling, c))[iv]
#        prob <- multCurves(dose, c)

#        return( -sum((resp2*weights2)*log(prob/(1-prob))+weights2*log(1-prob)) )
        return( -sum((resp*weights)[iv]*log(prob/(1-prob))+weights[iv]*log(1-prob)) )
    }    

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        total <- (object$"data")[iv, 5]
        success <- total*(object$"data")[iv, 2]    
        
        c(sum(log(choose(total, success))) - object$"fit"$"ovalue",  # object$"fit"$"ofvalue", 
        object$"sumList"$"lenData" - df.residual(object))  # object$"sumList"$"df.residual")
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }

#
#    ## Modifying ANOVA test (removing dose=0 and dose=Inf)
#    anovaTest2 <- function(formula, ds) {anovaTest(formula, ds[iv, ])}


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, 
    vcovfct = vcovfct, parmfct = parmfct))  # , anovaTest2=anovaTest2))
}


"drmLOFbinomial" <- function()
{
    ## Defining a goodness-of-fit test
    gofTest <- function(resp, weights, fitted, dfres)
    {
        ## Removing 0s and 1s in fitted values
        zeroTol <- 1e-12  # no global constant
        indVec <- ( (fitted < zeroTol) | (fitted > 1-zeroTol) )
        dfReduc <- sum(indVec)
    
        total <- weights  # (object$"data")[, 5]
        success <- resp*weights  # total*(object$"data")[, 2]
        expected <- total*fitted  # fitted(object) 
        
        ## Pearson's statistic (sum of squared Pearson residuals)
        c( sum( ((success - expected)^2 / (expected*(1 - fitted)))[!indVec] ), dfres - dfReduc)  # df.residual(object))    
    }


    ## Defining goodness-of-fit function
    anovaTest <- function(formula, ds)
    {
#       count <- resp*weights    
        anovaFit <- glm(formula, family=binomial(link = "logit"), data=ds)
        if (df.residual(anovaFit)>0) 
        {
            return(list(test = "lr", anovaFit = anovaFit))
        } else {
            return(NULL)
        }
    }
    anovaTest <- NULL  # lack-of-fit test not meaningful in most situations
    
    return(list(anovaTest = anovaTest, gofTest = gofTest))
}
