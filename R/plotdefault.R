"plotdefault" <- function(dose, resp, curve, dosePts, plotMat, colour=FALSE, conLevel=1e-2, conName="0", conNameYes, legend=TRUE, obs="average", 
                          col, log="x", lty, pch, xlab, ylab, xlim, ylim, ...) 
{

    ## Constructing the plot data

    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    }

    assayNo <- colConvert(curve)
    assayNoOld <- curve
    numAss <- length(unique(assayNo))

    numCol<-ncol(plotMat)

    maxR <- max(resp)
    options(warn=-1)  # suppressing warning in case maximum of NULL is taken
    maxPM <- apply(plotMat,2,max,na.rm=TRUE)
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    options(warn=0)

    xLimits <- xlim
    if (missing(ylim)) {yLimits <- c(min(resp),maxPM)} else {yLimits <- ylim}
    if (!missing(col)) {colourVec <- col}
    if (!missing(lty)) {ltyVec <- lty}    
    if (!missing(pch)) {pchVec <- pch}
#    print(pchVec)


#    ## Constructing vector of colours
#    colourVec <- rep(1, numCol)
#    if (is.logical(colour) && colour) {colourVec <- 1:numCol}
#    if (!is.logical(colour) && is.numeric(colour) && length(colour)==numCol) {colourVec <- colour}# else {stop("Argument 'colour' not correct")}
#    if (!missing(col)) {colourVec <- col}


    ## Plotting the fitted dose response curves
    par(las=1)
    
    if (!obs=="add")
    {    
    for (i in 1:numCol)
    {
        plotPoints <- switch(obs, "average" = cbind(as.numeric(names(tapVec<-tapply(resp[assayNo==i], dose[assayNo==i], mean))), tapVec),
                                  "none"    = c(max(dosePts)+10,max(c(maxPM,max(resp)))+10),
                                  "points"  = cbind(dose[assayNo==i], resp[assayNo==i]))
        
        if (i==1) 
        {        
            plot(plotPoints, xlab=xlab, ylab=ylab, log=log, xlim=xLimits, ylim=yLimits, pch=pchVec[i], axes=FALSE, frame.plot=TRUE, col=colourVec[i], ...) 
            
            xaxisTicks <- axTicks(1)             
            yaxisTicks <- axTicks(2)
            axis(2, at=yaxisTicks)
            xLabels <- as.character(xaxisTicks)
            if (conNameYes) {xLabels[1] <- conName}
            axis(1, at=xaxisTicks, labels=xLabels)

        } else {
            points(plotPoints, pch=pchVec[i], col=colourVec[i])
        }
    }
    }
    
    noPlot <- rep(FALSE, numCol)
    for (i in 1:numCol)
    {
        if (any(is.na(plotMat[,i]))) {noPlot[i] <- TRUE; next}
        lines(dosePts, plotMat[,i], lty=ltyVec[i], col=colourVec[i])  # solid only used if there are more curves than types of dashed lines
    }
    ltyIndex <- (1+1:numAss)
    ltyIndex[noPlot] <- 0
    
    if (legend && !(obs=="add"))
    {
        legend(xLimits[2], yLimits[2], as.character(unique(assayNoOld)), lty=ltyIndex, pch=1:numAss, col=colourVec[i], bty="n", xjust=1, yjust=1)
    }
    par(las=0)
    
    retData <- data.frame(dosePts, as.data.frame(plotMat))
    colnames(retData) <- c(xlab, as.character(unique(assayNoOld)))
    
    invisible(retData)
    
}
