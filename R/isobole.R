"isobole" <- 
function(object1, object2, exchange = 1, cifactor = 2, ename = "e", xaxis = "100",  
xlab, ylab, xlim, ylim, ...)
{

    parmVec <- coef(object1)
    namesPV <- names(parmVec)
    lenNPV <- length(namesPV)

    indVec <- regexpr(paste(ename, ":", sep = ""), namesPV, fixed = TRUE) > 0
    eVec <- parmVec[indVec]
    seVec <- (summary(object1)$"coefficients")[indVec, 2]

    edMat <- ED(object1, 50, display = FALSE)
    eVec <- (as.vector(edMat[, 1]))  # stripping off names
    seVec <- (as.vector(edMat[, 2]))  # stripping off names

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

#    if (identical(xaxis, "0")) 
#    {
#        eVec <- rev(eVec)
#        seVec <- rev(seVec)
#        mixProp <- rev(mixProp)
#    }

    Ex <- eVec * mixProp
    Ey <- eVec * (1-mixProp) * exchange

    lowerE <- eVec - cifactor*seVec
#    lowerEx <- lowerE*mixProp
#    lowerEx <- lowerE*cos(mixProp*pi/2)
#    lowerEy <- lowerE*(1-mixProp)*exchange    
#    lowerEy <- lowerE*sin(mixProp*pi/2)*exchange    

    upperE <- eVec + cifactor*seVec
#    upperEx <- upperE*mixProp
#    upperEx <- upperE*cos(mixProp*pi/2)
#    upperEy <- upperE*(1-mixProp)*exchange
#    upperEy <- upperE*sin(mixProp*pi/2)*exchange

#    lowerE <- eVec - 2 * seVec
    lowerEx <- lowerE * mixProp
    lowerEy <- lowerE * (1 - mixProp) * exchange
#    upperE <- eVec + 2 * seVec
    upperEx <- upperE * mixProp
    upperEy <- upperE * (1 - mixProp) * exchange


    ## Defining plotting frame
    if (missing(xlim)) {xLimits <- c(0, max(upperEx))} else {xLimits <- xlim}       
    if (missing(ylim)) {yLimits <- c(0, max(upperEy))} else {yLimits <- ylim}
    
    if (missing(xlab)) {xLab <- ifelse(identical(xaxis, "100"), "100", "0")} else {xLab <- xlab}
    if (missing(ylab)) {yLab <- ifelse(identical(xaxis, "100"), "0", "100")} else {yLab <- ylab}

    ## Swapping axes
    if (identical(xaxis, "0")) 
    {
        ETemp <- Ex
        Ex <- Ey
        Ey <- ETemp
        
        lowerTemp <- lowerEx
        lowerEx <- lowerEy
        lowerEy <- lowerTemp

        upperTemp <- upperEx
        upperEx <- upperEy
        upperEy <- upperTemp

        limTemp <- xLimits
        xLimits <- yLimits
        yLimits <- limTemp 
        
        mixProp <- 1 - mixProp
    }
    
    ## Drawing the plot frame    
    plot(0, type = "n", xlab = xLab, ylab = yLab, xlim = xLimits, ylim = yLimits, ...)  # empty plot
    
    ## Plotting rays in first quadrant
    raySlopes <- (1 - mixProp) / mixProp
#    raySlopes <- mixProp/(1 - mixProp) 
#    raySlopes <- tan(mixProp * pi/2)
    for (i in raySlopes[is.finite(raySlopes)]) {abline(0, exchange*i, lty = 3)}
    abline(v = 0, , lty = 3)  # adding vertical line (for infinite slope)       

#    raySlopes <- mixProp/(1 - mixProp)
#    for (i in raySlopes[is.finite(raySlopes)]) 
#    {
#        abline(0, exchange * i, lty = 3)
#    }

    ## Plotting ED50 values with confidence intervals  
    points(Ex, Ey, pch = 19)
    segments(lowerEx, lowerEy, upperEx, upperEy, lwd = 2)    
    
#    points(eVec*cos(mixProp*pi/2), eVec*sin(mixProp*pi/2)*exchange, pch = 19)

#    katx <- function(eVal, slope) {cos(atan(slope))*eVal}
#    katy <- function(eVal, slope) {sin(atan(slope))*eVal}    

#    points(katx(eVec, raySlopes), katy(eVec, raySlopes)*exchange, pch = 19)
#    segments(lowerEx, lowerEy, upperEx, upperEy, lwd = 2)
# old    segments(katx(lowerE, raySlopes), katy(lowerE, raySlopes)*exchange, 
#        katx(upperE, raySlopes), katy(upperE, raySlopes)*exchange, lwd = 2)

    if (!missing(object2))
    {

        ## Retrieving parameter estimates from fit of Hewlett's model
        parmVec <- coef(object2)
        namesPV <- names(parmVec)
        lenNPV <- length(namesPV)

        curveStr1 <- paste("I(1/(", object1$curveVarNam, "/100))", sep = "")
        curveStr2 <- paste("I(1/(1 - ", object1$curveVarNam, "/100))", sep = "")
        ED50.1 <- parmVec[regexpr(curveStr1, namesPV, fixed = TRUE) > 0]
        ED50.2 <- parmVec[regexpr(curveStr2, namesPV, fixed = TRUE) > 0]
        ## maybe use an extractor like 'EDmix' instead of the above 4 lines
        
#        ED50.1 <- parmVec[regexpr("I(1/(pct1/100))", namesPV, fixed = TRUE) > 0]
#        ED50.2 <- parmVec[regexpr("I(1/(1 - pct1/100))", namesPV, fixed = TRUE) > 0]
    
        ## Plotting the isobole based on the Voelund model
        ## Based on R lines from Helle Sørensen
        if (inherits(object2, "Voelund"))
        {
            Voelund <- function(t0, m0, eta, andeli)
            {
                if (andeli == 1) 
                {
                    tmp <- log(t0) 
                } else if (andeli == 0) 
                {
                    tmp <- log(m0)
                } else {
                    propor <- (1- andeli) / andeli
                    tmp <- (1+t0/m0*propor)^(1-eta[1]) + (t0/m0*propor)^eta[2] * (1+t0/m0*propor)^(1-eta[2])
                    tmp <- t0/tmp    
                    tmp <- log(tmp*(1+propor))
                } 
                tmp
            }
            Voelund2 <- Vectorize(Voelund, "andeli")
                                     
            mixVec <- seq(1, 0, length = 100)  # oops hardcoded "100"
            voeVec <- exp(Voelund2(ED50.1, ED50.2, tail(parmVec, 2), mixVec))        
            tvals <- mixVec * voeVec
            mvals <- (1 - mixVec) * voeVec
#            tvals <- voeVec * cos(mixVec*pi/2)
#            mvals <- voeVec * sin(mixVec*pi/2)

            xVal <- tvals
#            yVal <- tvals
            yVal <- exchange * mvals
#            xVal <- exchange * mvals

        } else {
    
            ## Plotting the isobole based on the Hewlett model
            ## Based on R lines from Helle Sørensen
            Hewlett <- function(t, t0, m0, lambda)
            {
                tmp <- (m0^(1/lambda) - (t*m0/t0)^(1/lambda))^lambda
                tmp[ t > t0 ] <- 0
                tmp
            }    
            if (inherits(object2, "Hewlett")) 
            {
                lambda <- parmVec[lenNPV]            
            } else {
                lambda <- 1
            }
#            yVal <- seq(0, ED50.1, length = 100)
#            xVal <- exchange*Hewlett(yVal, ED50.1, ED50.2, lambda)

            xVal <- seq(0, ED50.1, length = 100)
            yVal <- exchange * Hewlett(xVal, ED50.1, ED50.2, lambda)
        }
        if (identical(xaxis, "0")) 
        {
            tempVal <- xVal
            xVal <- yVal
            yVal <- tempVal
        }    
        lines(xVal, yVal)        
    }

#    invisible()
}
