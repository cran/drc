"mr.test" <- function(object1, object2, object, x, var.equal = TRUE, component = 1)
{
    dnFct <- deriv(~c+(d-c)/(1+(x/e)^b),c("b","c","d","e"), function(b,c,d,e,x){}, hessian = TRUE)
#    dnFct <- deriv(~d*exp(-e*x),c("d","e"), function(d,e,x){}, hessian = TRUE)  # exponential model
#    dnFct <- deriv(~d*exp(-e/x),c("d","e"), function(d,e,x){}, hessian = TRUE)  # exponential model
    daFct <- deriv(~c+(d-c)*exp(-exp(b*(log(x)-log(e)))),c("b","c","d","e"), function(b,c,d,e,x){}, 
    hessian = TRUE)

    beta1 <- coef(object1)
    beta2 <- coef(object2)
    beta <- coef(object)
    diffVec <- beta2 - beta
    
    lenx <- length(x)    
#    if (identical(var.equal, FALSE))
#    {
#        res1 <- residuals(object1)
#        res2 <- residuals(object2)
#    } else {
#        res1 <- rep(sqrt(summary(object1)$resVar), lenx) 
#        res2 <- res1 
#    }
#    res1 <- rep(1, lenx)
#    res2 <- res1
    
    Df1 <- attr(dnFct(beta1[1], beta1[2], beta1[3], beta1[4], x), "gradient")  # * res1
#    Df1 <- attr(dnFct(beta1[1], beta1[2], x), "gradient")  # exponential model
    
#    D2g1 <- attr(dnFct(theta1[1],theta1[2],theta1[3],theta1[4], x), "hessian")  # not needed
    Df1Df1 <- matrix(apply(apply(Df1, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
#    print(Df1)
#    print(Df1Df1)

    Df2 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "gradient")  # * res1
    Df2Df2 <- matrix(apply(apply(Df2, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
#    print(Df2)    
#    print(Df2Df2)

    ## Product of first derivatives from the two models
#    lenx <- length(x)
    Df1Df2.0 <- matrix(0, 4 * 4, lenx)
    for (i in 1:lenx)
    {
        Df1Df2.0[, i] <- as.vector(Df1[i, ]%*%t(Df2[i, ]))   
    }
    Df1Df2 <- matrix(apply(Df1Df2.0, 1, mean), 4, 4)

    ## Second derivatives for the alternative model
    fitDiff <- fitted(object2) - fitted(object1)
#    print(fitDiff)
    D2f2.0 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "hessian") 
    D2f2 <- matrix(0, 4, 4)
    for (i in 1:4)
    {
        D2f2[i, ] <- apply(D2f2.0[, , i] * fitDiff, 2, mean) 
    }
#    print(Df2Df2)
#    print(D2f2)

#    ## Using derivatives without residuals multiplied
#    Df2.2 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "gradient")
#    Df2Df2.2 <- matrix(apply(apply(Df2.2, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
    H <- Df2Df2 + D2f2
    Hinv <- solve(H)
    cov12 <- Df1Df2
    cov21 <- t(Df1Df2)
    var1 <- Df1Df1
    var2 <- Df2Df2

    Wfit <- Hinv %*% cov21 %*% solve(var1) %*% cov12 %*% Hinv
#    Wfit <- Hinv %*% cov21 %*% ginv(var1) %*% cov12 %*% Hinv  # exponential model
    Wobs <- Hinv %*% var2 %*% Hinv
    
#    varEst <- (Wobs - Wfit)/(lenx*summary(object1)$resVar)
    varEst <- (Wobs - Wfit) * (summary(object1)$resVar) / lenx
#    varEst <- (Wobs - Wfit)  / lenx

#    print(diffVec/sqrt(diag(varEst)))
#    print(varEst)
#    print(solve(varEst[1:4,1:4]))
     if (identical(var.equal, FALSE))
     {
         Hinv <- Hinv / lenx
         varEst0 <- - Hinv %*% t(Df2) + Hinv %*% cov21 %*% solve(var1) %*% t(Df1)
         
         p2 <- length(beta2)
         veMat <- matrix(NA, lenx, p2^2)
         res1 <- residuals(object1)
         for (i in 1:lenx)
         {
#             veMat[i, ] <- as.vector(outer(varEst0[, i], varEst0[, i])) * summary(object1)$resVar * lenx  # (res1[i]^2)
             veMat[i, ] <- as.vector(outer(varEst0[, i], varEst0[, i])) * (res1[i]^2) * lenx
         }
         varEst <- matrix(apply(veMat, 2, mean), p2, p2)
     }
    
#    sInd <- c(1, 2, 3, 4)
#    chiSq <- diffVec[sInd]%*%solve(varEst[sInd, sInd])%*%diffVec[sInd]
#    chiSq <- diffVec[sInd]%*%ginv(varEst[sInd, sInd])%*%diffVec[sInd]
#    pVal <- 1 - pchisq(chiSq, length(sInd))
    
    
    stdErr <- sqrt(varEst[component, component])
    diffVal <- diffVec[component]
    chiSq <- diffVal / stdErr
    pVal <- 2 * (1 - pnorm(abs(chiSq)))

    retVec <- c(chiSq, pVal, diffVal ,stdErr)
    names(retVec) <- c("Statistic", "p-value", "Difference", "SE")
    return(retVec)
}


"sim.mr.test" <- function(noSim, noObs = 20, noRep = NA, seed = 20070327, 
true = c("llogistic", "weibull", "expdecay"), var.equal = TRUE, 
xVec, sdVec, etmotc = FALSE, ...)
{
    true <- match.arg(true)
    
    set.seed(seed)

    expdecay <-
    function(x, b = 4, d = 2.15)
    {
        d * exp(-b * x)
    }    
    llogist <- 
    function(x, c = 0.08, d = 2.15, b = 1.84, e = 0.2)
    {
        c+(d-c)/(1+(x/e)^b)
    }
    weibull <-
    function(x, c = 0.08, d = 2.15, b = 1, e = 0.2)
    {
        c+(d-c)*exp(-exp(b*(log(x)-log(e))))
    }
    
    genFct <- switch(true, "expdecay" = expdecay, "llogistic" = llogist, "weibull" = weibull) 
    
    if (!is.na(noRep))
    {
        if (missing(xVec)) 
        {
            xFct <- function(){rep( (1:noObs)/(noObs+1), rep(noRep, noObs))}
        } else {
            xFct <- function(){xVec}
        }
#        errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054))}
        
        if (!missing(sdVec))
        {  # variance heterogeneity
#            errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054) * sdVec)}        
            errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054) * rep(3*exp(-c(1:10)/3), rep(noRep, noObs)))} # noObs=10           
        } else {  # variance homogeneity
            errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054))}
        }
    } else {
        if (missing(xVec)) 
        {
            xFct <- function(){sort(runif(noObs, 0, 1))}
        } else {
            xFct <- function(){xVec}
        }

        if (!missing(sdVec))
        {  # variance heterogeneity
            errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054) * sdVec)}                    
        } else {  # variance homogeneity
            errFct <- function(){rnorm(noObs, 0, sqrt(0.0054))}
        }
#        errFct <- function(){rnorm(noObs, 0, sqrt(0.0054))}
#        errFct <- function(){rnorm(noObs, 0, sqrt(0.0054) * (2*exp(-c(1:100)/70)) )}  # noObs = 100
    }

    if (etmotc)
    {
        errFct <- function(){rnorm(15, 0,  0.00737)}    
        xFct <- function(){c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        0.01, 0.36072, 0.90180, 2.70540, 5.41080, 14.42880, 36.07200, 90.18000)}
    
        genFct <-
        function(x, c = 0.22, d = 0.65, b = 0.94, e = 1.17)
        {
            c+(d-c)/(1+(x/e)^b)
        }    
    }   
    
    pVec <- rep(NA, noSim)
    seVec <- rep(NA, noSim)
    tsVec <- rep(NA, noSim)
    for (i in 1:noSim)
    {
#        x <- sort(runif(noObs, 0, 1))
        x <- xFct()
        y <- genFct(x) + errFct() 
#        print(cbind(x,y))

        if (!etmotc)
        {   
            fit1 <- try(multdrc(y ~ x, start = c(1.84, 0.08, 2.15, 0.2)), silent = TRUE)  # fitting log-logistic model
            fit2 <- try(multdrc(y ~ x, fct=w4(), start = c(1, 0.08, 2.15, 0.2)), silent = TRUE)  # fitting Weibull model
            fit <- try(multdrc(fitted(fit1) ~ x, fct=w4(), start = c(1, 0.08, 2.15, 0.2)), silent = TRUE)
        } else {
            fit1 <- try(multdrc(y ~ x, start = c(0.94, 0.22, 0.65, 1.17)), silent = TRUE)  # fitting log-logistic model
            fit2 <- try(multdrc(y ~ x, fct=w4(), start = c(-0.94, 0.22, 0.65, 1.17)), silent = TRUE)  # fitting Weibull model
            fit <- try(multdrc(fitted(fit1) ~ x, fct=w4(), start = c(-0.94, 0.22, 0.65, 1.17)), silent = TRUE)
        }
    
        if (inherits(fit1, "try-error") || inherits(fit2, "try-error") || inherits(fit, "try-error"))
        {
            pVec[i] <- NA
        } else {
            pVal <- try(mr.test(fit1, fit2, fit, x, var.equal = var.equal)[2], silent = TRUE)
            if (inherits(pVal, "try-error"))
            {
                pVec[i] <- NA
            } else {
                pVec[i] <- pVal
            }
        }

#        testRes <- mr.test(fit1, fit2, fit, x)        
#        pVec[i] <- testRes[2]
#        seVec[i] <- testRes[3]
#        tsVec[i] <- testRes[1]        
    }
    pVec
#    list(pVec, seVec, tsVec)
}

if (FALSE)
{
    ### Design without replicates  
    ## Examining the power of the test
    pmMR <- matrix(NA, 1000, 7)
    colnames(pmMR) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

    pmMR[, 1] <- sim.mr.test(1000, 10, seed=2221, true = "weibull")
    pmMR[, 2] <- sim.mr.test(1000, 25, seed=2222, true = "weibull")
    pmMR[, 3] <- sim.mr.test(1000, 50, seed=2223, true = "weibull")
    pmMR[, 4] <- sim.mr.test(1000, 75, seed=2224, true = "weibull")
    pmMR[, 5] <- sim.mr.test(1000, 100, seed=2225, true = "weibull")
    pmMR[, 6] <- sim.mr.test(1000, 250, seed=2226, true = "weibull")
    pmMR[, 7] <- sim.mr.test(1000, 500, seed=2227, true = "weibull")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmMR, 2, pFct)

    ## Examining the power of the test
    ##  under missspecified alternative
    pmMR.2 <- matrix(NA, 1000, 7)
    colnames(pmMR) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

    pmMR.2[, 1] <- sim.mr.test(1000, 10, seed=2221, true = "expdecay")
    pmMR.2[, 2] <- sim.mr.test(1000, 25, seed=2222, true = "expdecay")
    pmMR.2[, 3] <- sim.mr.test(1000, 50, seed=2223, true = "expdecay")
    pmMR.2[, 4] <- sim.mr.test(1000, 75, seed=2224, true = "expdecay")
    pmMR.2[, 5] <- sim.mr.test(1000, 100, seed=2225, true = "expdecay")
    pmMR.2[, 6] <- sim.mr.test(1000, 250, seed=2226, true = "expdecay")
#    pmMR[, 7] <- sim.mr.test(1000, 500, seed=2227, true = "weibull")
    pFct2 <- function(x, alpha = 0.05) {sum(x < alpha, na.rm = TRUE)/1000}    
    apply(pmMR.2, 2, pFct2)

    ## Examining the level of the test
    pmMRt <- matrix(NA, 1000, 7)
    colnames(pmMRt) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
    
    pmMRt[, 1] <- sim.mr.test(1000, 10, seed = 2231)
    pmMRt[, 2] <- sim.mr.test(1000, 25, seed = 2232)
    pmMRt[, 3] <- sim.mr.test(1000, 50, seed = 2233)
    pmMRt[, 4] <- sim.mr.test(1000, 75, seed = 2234)
    pmMRt[, 5] <- sim.mr.test(1000, 100, seed = 2235)
    pmMRt[, 6] <- sim.mr.test(1000, 250, seed = 2236)
    pmMRt[, 7] <- sim.mr.test(1000, 500, seed = 2237)
    apply(pmMRt, 2, pFct)

    ## Without assuming variance homogeneity
    pmMRt.2 <- matrix(NA, 1000, 7)
    colnames(pmMRt.2) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
    
    pmMRt.2[, 1] <- sim.mr.test(1000, 10, seed = 2231, var.equal = FALSE)
    pmMRt.2[, 2] <- sim.mr.test(1000, 25, seed = 2232, var.equal = FALSE)
    pmMRt.2[, 3] <- sim.mr.test(1000, 50, seed = 2233, var.equal = FALSE)
    pmMRt.2[, 4] <- sim.mr.test(1000, 75, seed = 2234, var.equal = FALSE)
    pmMRt.2[, 5] <- sim.mr.test(1000, 100, seed = 2235, var.equal = FALSE)
    pmMRt.2[, 6] <- sim.mr.test(1000, 250, seed = 2236, var.equal = FALSE)
    pmMRt.2[, 7] <- sim.mr.test(1000, 500, seed = 2237, var.equal = FALSE)
    apply(pmMRt.2, 2, pFct)



    ### Design with replicates
    ## Examining the power of the test
    pmMR.rep <- matrix(NA, 1000, 6)
    colnames(pmMR.rep) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmMR.rep[, 1] <- sim.mr.test(1000, 10, noRep = 1, seed=200805071, true = "weibull")
    pmMR.rep[, 2] <- sim.mr.test(1000, 10, noRep = 2, seed=200805072, true = "weibull")
    pmMR.rep[, 3] <- sim.mr.test(1000, 10, noRep = 3, seed=200805073, true = "weibull")
    pmMR.rep[, 4] <- sim.mr.test(1000, 10, noRep = 5, seed=200805074, true = "weibull")
    pmMR.rep[, 5] <- sim.mr.test(1000, 10, noRep = 10, seed=200805075, true = "weibull")
    pmMR.rep[, 6] <- sim.mr.test(1000, 10, noRep = 20, seed=200805076, true = "weibull")
    apply(pmMR.rep, 2, pFct)

    ## Examining the power of the test
    ##  under misspecified alternative
    pmMR.rep.2 <- matrix(NA, 1000, 6)
    colnames(pmMR.rep.2) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmMR.rep.2[, 1] <- sim.mr.test(1000, 10, noRep = 1, seed=200805141, true = "expdecay")
    pmMR.rep.2[, 2] <- sim.mr.test(1000, 10, noRep = 2, seed=200805142, true = "expdecay")
    pmMR.rep.2[, 3] <- sim.mr.test(1000, 10, noRep = 3, seed=200805143, true = "expdecay")
    pmMR.rep.2[, 4] <- sim.mr.test(1000, 10, noRep = 5, seed=200805144, true = "expdecay")
    pmMR.rep.2[, 5] <- sim.mr.test(1000, 10, noRep = 10, seed=200805145, true = "expdecay")
    pmMR.rep.2[, 6] <- sim.mr.test(1000, 10, noRep = 20, seed=200805146, true = "expdecay")
    apply(pmMR.rep.2, 2, pFct2)

    ## Under variance heterogeneity
    pmMR.rep.3 <- matrix(NA, 1000, 5)
    colnames(pmMR.rep.3) <- as.character(c(2, 3, 5, 10, 20))

    pmMR.rep.3[, 1] <- sim.mr.test(1000, 10, noRep = 2, seed=200805162, true = "weibull", sdVec = NA, var.equal = FALSE)
    pmMR.rep.3[, 2] <- sim.mr.test(1000, 10, noRep = 3, seed=200805163, true = "weibull", sdVec = NA, var.equal = FALSE)
    pmMR.rep.3[, 3] <- sim.mr.test(1000, 10, noRep = 5, seed=200805164, true = "weibull", sdVec = NA, var.equal = FALSE)
    pmMR.rep.3[, 4] <- sim.mr.test(1000, 10, noRep = 10, seed=200805165, true = "weibull", sdVec = NA, var.equal = FALSE)
    pmMR.rep.3[, 5] <- sim.mr.test(1000, 10, noRep = 20, seed=200805166, true = "weibull", sdVec = NA, var.equal = FALSE)
    apply(pmMR.rep.3, 2, pFct2)

    ## Examining the level of the test
    pmMRt.rep <- matrix(NA, 1000, 6)
    colnames(pmMRt.rep) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmMRt.rep[, 1] <- sim.mr.test(1000, 10, noRep = 1, seed=2231)    
#    pmMRt.rep[, 1] <- sim.mr.test(1000, 10, noRep = 1, seed=200805081)
#    pmMRt.rep[, 2] <- sim.mr.test(1000, 10, noRep = 2, seed=200805082)
    pmMRt.rep[, 2] <- sim.mr.test(1000, 10, noRep = 2, seed=200805133)
    pmMRt.rep[, 3] <- sim.mr.test(1000, 10, noRep = 3, seed=2008050813)
    pmMRt.rep[, 4] <- sim.mr.test(1000, 10, noRep = 5, seed=2008050814)
    pmMRt.rep[, 5] <- sim.mr.test(1000, 10, noRep = 10, seed=2008050815)
    pmMRt.rep[, 6] <- sim.mr.test(1000, 10, noRep = 20, seed=2008050816)
    apply(pmMRt.rep, 2, pFct)

    ## Under variance heterogeneity
    pmMRt.rep.3 <- matrix(NA, 1000, 5)
    colnames(pmMRt.rep.3) <- as.character(c(2, 3, 5, 10, 20))
    pmMRt.rep.3[, 1] <- sim.mr.test(1000, 10, noRep = 2, seed=200805162, true = "llogistic", sdVec = NA, var.equal = FALSE)
    pmMRt.rep.3[, 2] <- sim.mr.test(1000, 10, noRep = 3, seed=200805163, true = "llogistic", sdVec = NA, var.equal = FALSE)
    pmMRt.rep.3[, 3] <- sim.mr.test(1000, 10, noRep = 5, seed=200805164, true = "llogistic", sdVec = NA, var.equal = FALSE)
    pmMRt.rep.3[, 4] <- sim.mr.test(1000, 10, noRep = 10, seed=200805165, true = "llogistic", sdVec = NA, var.equal = FALSE)
    pmMRt.rep.3[, 5] <- sim.mr.test(1000, 10, noRep = 20, seed=200805166, true = "llogistic", sdVec = NA, var.equal = FALSE)
    apply(pmMRt.rep.3, 2, pFct2)

        
}
