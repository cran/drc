"boxcox.drc" <- function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
eps = 1/50, bcAdd = 0, level = 0.95, xlab = expression(lambda), ylab = "log-Likelihood", ...)
{

    lenlam <- length(lambda)
    llVec <- rep(NA, lenlam)
    for (i in 1:lenlam)
    {
        drcTemp <- try(update(object, bc = lambda[i], bcAdd = bcAdd), silent = TRUE)
        if (!inherits(drcTemp, "try-error")) {llVec[i] <- logLik(drcTemp)}
    }
    lv <- lambda[which.max(llVec)]
    llv <- max(llVec, na.rm = TRUE)
    ci <- boxcoxCI(lambda, llVec, level)

    if (plotit)  # based on boxcox.default
    {
        plot(lambda, llVec, type ="l", xlab = xlab, ylab = ylab, ...)
        
        plims <- par("usr")
        y0 <- plims[3]
        lim <- llv-qchisq(level, 1)/2
        
        segments(lv, llv, lv, y0, lty=3)
        segments(ci[1], lim, ci[1], y0, lty = 3)  # lower limit
        segments(ci[2], lim, ci[2], y0, lty = 3)  # upper limit
        
        scal <- (1/10 * (plims[4] - y0))/par("pin")[2] 
        scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1] 
        text(lambda[1] + scx, lim + scal, " 95%") 
        abline(h = lim, lty = 3)

    } 
    retFit <- update(object, bc = lv, bcAdd = bcAdd)
    retFit$lambda <- list(lambda = lv, ci = ci)
    retFit$boxcox[c(2, 3)] <- ci
    ## future: make boxcox and lambda into one component in the fit
     
    invisible(retFit)
}

"boxcoxCI" <- 
function(x, y, level = 0.95)
{
    ## R lines taken from boxcox.default in the package MASS and then slightly modified
    xl <- x
    loglik <- y
    
    llnotna <- !is.na(loglik)
    xl <- xl[llnotna]
    loglik <- loglik[llnotna]
    
    m <- length(loglik)

    mx <- (1:m)[loglik == max(loglik)][1]
    Lmax <- loglik[mx]
    lim <- Lmax - qchisq(level, 1)/2

    ind <- range((1:m)[loglik > lim])

    xx <- rep(NA, 2)
    if(loglik[1] < lim) 
    {
        i <- ind[1]
        xx[1] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])

    }
    if(loglik[m] < lim) 
    {
        i <- ind[2] + 1
        xx[2] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
    }
    return(xx)
}


