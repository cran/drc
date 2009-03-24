"llogistic.ssf" <- function(method = c("1", "2", "3"), fixed)
{
    method <- match.arg(method)
    numParm <- length(fixed)

    notFixed <- is.na(fixed)

    ## Version 1 (default)    
    ssFct <- switch(method,
    "1" = 
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        yRange <- range(y)
        lenyRange <- 0.001 * diff(yRange)
        cVal <- yRange[1] - lenyRange  # the c parameter        
        dVal <- yRange[2] + lenyRange  # the d parameter
        fVal <- 1  # better choice may be possible!        

        beVal <- find.be2(x, y, cVal, dVal)
        
        return(c(beVal[1], cVal, dVal, beVal[2], fVal)[notFixed])
    },
    "2" =
    function(dframe, doseScaling, respScaling)
    {
        x <- dframe[, 1] / doseScaling
        y <- dframe[, 2] / respScaling

        cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2] / respScaling)
        dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3] / respScaling)

        fVal <- 1  # need not be updated with value in 'fixed[5]'
        # better choice than 1 may be possible! 
        # the f parameter, however, is very rarely a magnitude of 10 larger or smaller

#        # Cutting away response values close to d
#        indexT1a <- x > 0
#        x2 <- x[indexT1a]
#        y2 <- y[indexT1a]
            
        beVal <- find.be2(x, y, cVal, dVal)
# These lines are not needed as the b and e parameters are not used in further calculations
#        bVal <- ifelse(notFixed[1], beVec[1], fixed[1])
#        eVal <- ifelse(notFixed[4], beVec[2], fixed[4] / doseScaling)

#        bVal <- beVec[1]
#        eVal <- beVec[2]                
#        return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        
        return(c(beVal[1], cVal, dVal, beVal[2], fVal)[notFixed])
    },
    "3" =
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2])
        dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3])
        fVal <- 1  # need not be updated with value in 'fixed[5]'
        
#        if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
#        # only estimate of upper limit if a single unique dose value 
# no longer needed
           
        beVec <- find.be1(x, y, cVal, dVal)
        bVal <- beVec[1]
        eVal <- beVec[2]
            
        return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
    })

    ## Finding b and e based on stepwise increments
    find.be1 <- function(x, y, cVal, dVal)
    {
        unix <- unique(x)
        uniy <- tapply(y, x, mean)
        lenx <- length(unix)
        
        j <- 2
        for (i in 2:lenx)
        {
            crit1 <- (uniy[i] > (cVal + dVal)/2) && (uniy[i-1] < (cVal + dVal)/2)
            crit2 <- (uniy[i] < (cVal + dVal)/2) && (uniy[i-1] > (cVal + dVal)/2)
            if (crit1 || crit2) break
            j <- j + 1
        }
        eVal <- (unix[j] + unix[j-1]) / 2
        bVal <- -sign(uniy[j] - uniy[j-1])  # -(uniy[j] - uniy[j-1]) / (unix[j] - unix[j-1])
        return(as.vector(c(bVal, eVal)))  
    }
    
#    ## Finding b and e based on linear regression
#    find.be2 <- function(x, y, c, d)
#    {
#        logitTrans <- log((d - y)/(y - c))  
#
#        lmFit <- lm(logitTrans ~ log(x))
##        eVal <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))
##        bVal <- coef(logitFit)[2]
#
#        coefVec <- coef(lmFit)
#        bVal <- coefVec[2]        
#        eVal <- exp(-coefVec[1]/bVal)    
#
#        return(as.vector(c(bVal, eVal)))
#    }
    
    ## Finding b and e based on linear regression
    find.be2 <- function(x, y, cVal, dVal)
    {
#        logitTrans <- log((d - y)/(y - c))  
        lmFit <- lm(log((dVal - y)/(y - cVal)) ~ log(x), subset = x > 0)
        coefVec <- coef(lmFit)
        bVal <- coefVec[2]        
        eVal <- exp(-coefVec[1] / bVal)

        return(as.vector(c(bVal, eVal)))
    }  

    ## Returning self starter function
    ssFct
}

LL.ssf <- llogistic.ssf