"isobole" <- 
function(object1, object2, ename = "e", xaxis = "0", exchange, xlab = "A", ylab = "B", xlim, ylim, ...)
{

    parmVec <- coef(object1)
    namesPV <- names(parmVec)
    lenNPV <- length(namesPV)

    indVec <- regexpr(paste(ename, ":", sep = ""), namesPV, fixed = TRUE) > 0
    eVec <- parmVec[indVec]
    seVec <- (summary(object1)$"coefficients")[indVec, 2]

    mixProp <- unique(object1$data[,4])
    mixProp <- mixProp[ (mixProp >= 0) & (mixProp <= 100) ]/100  # removing control level

#    posOnRay <- function(len, slope)
#    {
#        xPos <- sqrt( len^2 / (1+slope^2) )
#        yPos <- slope*xPos
#        
#        yPos[!is.finite(slope)] <- len[!is.finite(slope)]
#        
#        list(xPos, yPos)
#    }


    if (xaxis == "0") {eVec <- rev(eVec)}

    lowerE <- eVec - 2*seVec
    lowerEx <- lowerE*mixProp
    lowerEy <- lowerE*(1-mixProp)*exchange    
    upperE <- eVec + 2*seVec
    upperEx <- upperE*mixProp
    upperEy <- upperE*(1-mixProp)*exchange


    ## Defining plotting frame
    if (missing(xlim)) {xLimits <- c(0, max(upperEx))} else {xLimits <- xlim}       
    if (missing(ylim)) {yLimits <- c(0, max(upperEy))} else {yLimits <- ylim}    
    plot(0, type = "n", xlab = xlab, ylab = ylab, xlim = xLimits, ylim = yLimits, ...) 
    
    ## Plotting rays in first quadrant
    raySlopes <- mixProp/(1 - mixProp)   
    for (i in raySlopes[is.finite(raySlopes)]) {abline(0, exchange*i, lty = 3)}
    abline(v = 0, , lty = 3)  # adding vertical line (for infinite slope)       

    ## Plotting ED50 values with confidence intervals
    points(eVec*mixProp, eVec*(1-mixProp)*exchange, pch = 19)
    segments(lowerEx, lowerEy, upperEx, upperEy, lwd=2)


    if (!missing(object2))
    {

        ## Retrieving parameter estimates from fit of Hewlett's model
        parmVec <- coef(object2)
        namesPV <- names(parmVec)
        lenNPV <- length(namesPV)

        ED50.1 <- parmVec[regexpr("I(1/(pct1/100))", namesPV, fixed = TRUE) > 0]
        ED50.2 <- parmVec[regexpr("I(1/(1 - pct1/100))", namesPV, fixed = TRUE) > 0]
    
        if (inherits(object2, "Hewlett")) 
        {
            lambda <- parmVec[lenNPV]            
        } else {
            lambda <- 1
        }

        ## Hewlett-modellens isobol
        Hewlett <- function(t,t0,m0,lambda)
        {
            tmp <- (m0^(1/lambda) - (t*m0/t0)^(1/lambda))^lambda
            tmp[ t > t0 ] <- 0
            tmp
        }
            
        xVal <- seq(0, ED50.1, length = 100)
        yVal <- exchange*Hewlett(xVal, ED50.1, ED50.2, lambda)
            
        lines(xVal, yVal)
    }

#    invisible()
}
