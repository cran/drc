"plotbin" <- 
function(x, ..., colour=FALSE, conLevel=1e-2, conName="0", grid=100, legend=FALSE, obs="average", col, log="x", lty, pch, xlab, ylab, xlim, ylim)
{
    object <- x

    curveid <- object$data$curve
    dose <- object$data$dose
    resp <- object$data$resp  # proportions
    total <- object$data$total

    int <- object$data$int
    slope <- object$data$slope
    lower <- object$data$lower
    upper <- object$data$upper

    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    }

    ciNum <- colConvert(curveid)
    uniID <- unique(curveid)
    lenuid <- length(uniID)
    
    intLevels <- tapply(int, curveid, function(x) x[1])
    sloLevels <- tapply(slope, curveid, function(x) x[1])
    if (!is.null(lower)) {lowLevels <- tapply(lower, curveid, function(x) x[1])} else {lowLevels <- rep(NA, lenuid)}
    if (!is.null(upper)) {upLevels <- tapply(upper, curveid, function(x) x[1])} else {upLevels <- rep(NA, lenuid)}   
#    if (!is.null(fixedLow)) {flLevels <- tapply(fixedLow, curveid, function(x) x[1])}    
#    if (!is.null(fixedUp)) {fuLevels <- tapply(fixedUp, curveid, function(x) x[1])}    
#    print(intLevels)
#    print(sloLevels)
#    print(lowLevels)
#    print(upLevels)


    ## Determining range of dose values
    if (missing(xlim))
    {
        xLimits <- c(min(dose),max(dose))
    } else {
        xLimits <- xlim; # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
    }


    ## Handling small dose values
    conNameYes <- FALSE
    if (xLimits[1]<conLevel)
    {
        xLimits[1] <- conLevel
        smallDoses <- dose<conLevel
        dose[smallDoses] <- conLevel
        conNameYes <- TRUE
    }
    if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}


    ## Constructing dose values for plotting
    dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), length=grid))


    ## Finding maximum on response scale
    plotMat <- matrix(0, grid, lenuid)
    for (i in 1:lenuid)
    {
#        print(c(intLevels[i], sloLevels[i], lowLevels[i], upLevels[i]))
        plotMat[,i] <- predict(object, data.frame(dose = dosePts, int = rep(intLevels[i], grid), slope = rep(sloLevels[i], grid), 
              low = rep(lowLevels[i], grid), up = rep(upLevels[i], grid)))
    }


    ## Setting up plot parameters
    if (missing(log)) {logScale <- "x"} else {logScale <- log}
    if (missing(pch)) {pchVec <- unique(ciNum)} else {pchVec <- pch}
    nv <- object$names
    if (missing(xlab)) {xLab <- as.character(object$call[[2]][3])} else {xLab <- xlab}
    if (missing(ylab)) {yLab <- as.character(object$call[[2]][2])} else {yLab <- ylab}

    if (!colour) {colVec <- rep(1, lenuid)}
    if (colour)
    {
        if (missing(col))
        {
            colVec <- (1:lenuid) + 1  # +1 is to avoid a black curve
        } else {
            colVec <- col
        }
    }

    if (missing(lty))
    {
        ltyVec <- (1:lenuid) + 1  # +1 is to avoid a solid line (as first line drawn!)
    } else {
        ltyVec <- lty
    }


    ## Calling plotdefault
        
    plotdefault(dose, resp, curveid ,dosePts, plotMat, colour=colour, conLevel=conLevel, conName=conName, conNameYes=as.logical(conNameYes), legend=legend, 
                obs=obs, col=colVec, log=log, lty=ltyVec, pch=pchVec, xlab=xLab, ylab=yLab, xlim=xLimits, ylim, ...)     

}
