"drmEMWadley" <- 
function(dose, resp, multCurves, startVec, zeroTol = 1e-12, doseScaling = 1)
{

    ## Finding indices for doses that give contribution to likelihood function
#    iv <- ( (multCurves(dose, startVec) > zeroTol) & (multCurves(dose, startVec) < 1-zeroTol) )


    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {
        prob <- multCurves(dose / doseScaling, head(c, -1))
        omZT <- 1 - zeroTol
        prob[prob > omZT] <- omZT
        prob[prob < zeroTol] <- zeroTol
        nyVal <- tail(c, 1)
        nyVal * sum(1 - prob) - log(nyVal) * sum((resp)) - sum(resp * (1 - prob))        
    }    

    if (FALSE)
    {
    
    opfct<-function(c){
        lnyVal <- c[3]
#        probVal <- exp(c[1]+c[2]*log(doseVal))/(1+exp(c[1]+c[2]*log(doseVal)))
        probVal <- pnorm(c[1]+c[2]*log(doseVal))
        exp(lnyVal) * sum(1 - probVal) - lnyVal * sum((respVal)) - sum(respVal * (1 - probVal))
    }
    
    }
    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )
        
        c(-object$"fit"$value + sum(log(gamma(resp+1))))  # adding scale constant
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


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFWadley" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}


if (FALSE)
{
    wadley1949.1 <- data.frame(
    time = c(0.00000, 19.99862, 29.99163, 39.99447, 50.00345),
    resp = c(280, 233, 85, 29, 3)
    )

    wad.m1 <- drm(resp ~ time, data = wadley1949.1, fct = LN.3(), type="Poisson")
    summary(wad.m1)                                      
    log10(coef(wad.m1)[3])
    ## log ED50 is 1.417 in Finney)

    ## upper limit in Finney: 282.02 and the same in "drc"

    ## slope in Finney: 7.5034 (using log10)
    -7.5034 / log(10)
    ## so all estimates agree
   
    ## Calculating the deviance (2.77 in Finney)
    sum(residuals(wad.m1)^2 / (fitted(wad.m1)))


    S.capricoruntum <- data.frame(
    conc = c(0,0,0,0,0,1,1,1,1,1,5,5,5,5,5,10,10,10,10,10,50,50,50,50,50), 
    resp = c(219, 228, 202, 237, 228, 167, 158, 158, 175,167,105,123,105,105,105,88,88,61,61,88,61,44,35,35,44)
    )
    S.capricoruntum

    s.cap.m1 <- drm(resp~conc,  data=S.capricoruntum, fct=LL.3(), type="Poisson")
    summary(s.cap.m1)
    
    s.cap.m2 <- drm(resp~conc,  data=S.capricoruntum, fct=LN.3(), type="Poisson")
    summary(s.cap.m2)
    
    # Calculating the deviance
    sum(residuals(s.cap.m2)^2 / (fitted(s.cap.m2))) 
    # equal to 34.34 close to 34.21 reported by Morgan

    ## Reducing the relative tolerance to get closer to the results obtained by Morgan
    s.cap.m2b <- drm(resp~conc,  data=S.capricoruntum, fct=LN.3(), type="Poisson", control=drmc(relTol=1e-12))
    summary(s.cap.m2b)
    sum(residuals(s.cap.m2b)^2 / (fitted(s.cap.m2b)))


}