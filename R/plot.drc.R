"plot.drc" <-
function(x, ..., level = NULL, breakCurve = FALSE, col = FALSE, conLevel, conName, 
         grid = 100, legend = TRUE, legendText, type = "average", obs, lty, log = "x", pch, xlab, ylab, xlim, ylim)
{    
#    breakCurve <- FALSE
    colour <- col  # replace 'colour' throughout the file

    if (!missing(obs)) {stop("Use 'type' instead of 'obs'")}  # remove in the future

    object <- x
#    zeroEps <- 1e-4 # to avoid log(0) at lower limit of dose range in plot 

    if (inherits(object, "bindrc"))
    {
        plotbin(x, ..., colour = colour, conLevel = conLevel, conName = conName, grid = grid, legend = legend, obs = type, col = col, log = log, 
                xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
    } else {


    ## Constructing the plot data
    origData<-object$"data"
    doseDim <- ncol(origData) - 4  # subtracting 4 because the data frame also contains columns with response, 
                                   # curve no. (in new and old enumeration) and weights
    if (doseDim > 1) {stop("No plot features for plots in more than two dimensions")}

#    dose <- origData[, 1:doseDim] + zeroEps
    dose <- origData[, 1:doseDim]
    resp <- origData[, doseDim+1]

    assayNo <- origData[, 3]
    assayNoOld <- origData[, 4]
    numAss <- length(unique(assayNo))
    doPlot <- is.null(level) || any(unique(assayNoOld)%in%level)
    if (!doPlot) {stop("Nothing to plot")}

    plotFct <- (object$"curve")[[1]]
    logDose <- (object$"curve")[[2]]


    if (missing(conLevel)) 
    {
        conLevel <- ifelse(is.null(logDose), 1e-2, log(1e-2))
    }
    if (missing(conName)) 
    {
        if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
    }


    ## Assigning axis names
    varNames <- colnames(origData)[1:(doseDim+1)]  # axis names are the names of the dose variable and response variable in the original data set
    if (missing(xlab)) {if (varNames[1] == "") {xlab <- "Dose"} else {xlab <- varNames[1]}}
    if (missing(ylab)) {if (varNames[2] == "") {ylab <- "Response"} else {ylab <- varNames[2]}}
#    if (missing(log)) {log <- "x"}
      
      
    ## Determining range of dose values
    if (missing(xlim)) 
    {
        xLimits <- c(min(dose),max(dose))
    } else {
        xLimits <- xlim; # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
    }


    ## Handling small dose values
    conNameYes <- FALSE
    
#    if (logDose) {conLevel <- log(conLevel)}
    if (xLimits[1] < conLevel) 
    {
        xLimits[1] <- conLevel
        smallDoses <- dose<conLevel
        dose[smallDoses] <- conLevel
        conNameYes <- TRUE 
    }
    if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}


    ## Constructing dose values for plotting
    if (doseDim == 1) 
    {
        if (is.null(logDose))
        {
           dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), length=grid))
        } else {
           dosePts <- seq(xLimits[1], xLimits[2], length=grid)
        }
    } else {}  # No handling of multi-dimensional dose values


    ## Finding maximum on response scale
    if (is.null(logDose)) {plotMat <- plotFct(dosePts)} else {plotMat <- plotFct(logDose^(dosePts))}
    numCol <- ncol(plotMat)

    maxR <- max(resp)
    options(warn = -1)  # suppressing warning in case maximum of NULL is taken
    maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    options(warn=0)  

    if (missing(ylim)) {yLimits <- c(min(resp),maxPM)} else {yLimits <- ylim}


    ## Constructing vector of colours
    colourVec <- rep(1, numCol)
    if (is.logical(colour) && colour) {colourVec <- 1:numCol}  
    if (!is.logical(colour) && is.numeric(colour) && length(colour)==numCol) {colourVec <- colour}# else {stop("Argument 'colour' not correct")}
    if (!is.logical(colour) && is.numeric(colour) && (!(length(colour)==numCol)) ) {colourVec <- rep(colour, numCol)}
#    if (!missing(col)) {colourVec <- col}
 
      
    ## Defining line types  
#    ltyVec <- rep(1, numCol)
#    if (is.logical(colour) && colour) {colourVec <- 1:numCol}  
#    if (!is.logical(colour) && is.numeric(colour) && length(colour)==numCol) {colourVec <- colour}# else {stop("Argument 'colour' not correct")}
#    if (!is.logical(colour) && is.numeric(colour) && (!(length(colour)==numCol)) ) {colourVec <- rep(colour, numCol)}
    if (!missing(lty)) 
    {
        if (length(lty)==1) {ltyVec <- rep(lty, numCol)} else {ltyVec <- lty}
    } else {
        ltyVec <- 1:numCol
    }

    ## Defining plot characters 
    if (!missing(pch)) 
    {
        if (length(pch)==1) {pchVec <- rep(pch, numCol)} else {pchVec <- pch}
    } else {
        pchVec <- 1:numCol
    }


    ## Plotting the fitted dose response curves
    par(las=1)
    if (!is.null(logDose)) {if (log=="x") {log <- ""}; if ( (log=="xy") || (log=="yx") ) {log <- "y"}}
    
    if (!type=="add")
    {    
    for (i in 1:numCol)
    {
        plotPoints <- switch(type, "average" = cbind(as.numeric(names(tapVec<-tapply(resp[assayNo==i], dose[assayNo==i], mean))), tapVec),
                                   "none"    = c(max(dosePts)+10,max(c(maxPM,max(resp)))+10),
                                   "points"  = cbind(dose[assayNo==i], resp[assayNo==i]))
        
#        if (colour) {j <- i} else {j <- 1}  
        if (i==1) 
        {        
#            if (obs)
#            {
#                plot(dose[assayNo==i],resp[assayNo==i],xlab=xlab,ylab=ylab,log=log,xlim=xLimits,ylim=yLimits,pch=i,axes=FALSE,frame.plot=TRUE,...) 
#            } else {
#                plot(max(dosePts)+10,max(c(maxPM,max(resp)))+10,xlab=xlab,ylab=ylab,log=log,xlim=xLimits,ylim=yLimits,pch=i,axes=FALSE,frame.plot=TRUE,...) 
#                # an empty plot
#            }
     
            if (is.null(level) || unique(assayNoOld)[i]%in%level)
            {           
                 if (breakCurve)
                 {
#                print(xLimits)
                xLimits[1] <- xLimits[1]*0.5
#                print(xLimits)
#                print(plotPoints)

                plotPoints[plotPoints[, 1]==conLevel, 1] <- xLimits[1]
#                print(plotPoints)
                }

#                if (!type=="add")
#                {                    
                    plot(plotPoints, xlab = xlab, ylab = ylab, log = log, xlim = xLimits, ylim = yLimits, axes = FALSE, 
                         frame.plot = TRUE, col = colourVec[i], pch = pchVec[i], ...) 
                
                xaxisTicks <- axTicks(1)
                if (breakCurve) {xaxisTicks <- unique(c(xLimits[1], xaxisTicks))}
                             
                yaxisTicks <- axTicks(2)
                axis(2, at=yaxisTicks)
#                xLabels <- as.character(xaxisTicks)
                xLabels <- as.expression(xaxisTicks)
                if (conNameYes) {xLabels[1] <- conName}
                axis(1, at=xaxisTicks, labels=xLabels)
#                print(xaxisTicks)

#                 }
#                } else {
#                    points(plotPoints, pch=i, col=colourVec[i], ...)
#                }

                if (breakCurve)
                {
                    intDiff <- diff(log(xaxisTicks[1:2], 10))
                    
                    vps<-baseViewports()
                    pushViewport(vps$inner, vps$figure, vps$plot)
                    
                    xPoints <- log(10^(log(xaxisTicks[1], 10) + (c(1:5)/6)*intDiff))  # converting to natural logarithm
                    yPoints <- c(-0.8, 0.8)*(unit(yaxisTicks[1], "native") - unit(0, "npc"))
                    
                    print(xPoints)
                    print(yPoints) 
                     
                    grid.move.to(unit(xPoints[1], "native"), unit(-0.5, "npc"))
                    grid.line.to(unit(xPoints[3], "native"), unit(0.5, "npc"), gp=gpar(lwd=1, lty=1))

                    grid.move.to(unit(xPoints[3], "native"), unit(-0.5, "npc"))
                    grid.line.to(unit(xPoints[5], "native"), unit(0.5, "npc"), gp=gpar(lwd=1, lty=1))
                    
                    popViewport(0)
                }

            }
            

        } else {
#            if (obs)
#            {
#                points(dose[assayNo==i], resp[assayNo==i], pch=i, ...)

            if (breakCurve) 
            {
                xLimits[1] <- xLimits[1]*0.5            
                plotPoints[plotPoints[, 1]==conLevel, 1] <- xLimits[1]
            }

            matchLevel <- match(unique(assayNoOld)[i],level)
            if ( (!is.null(level)) && (!is.na(matchLevel)) && (matchLevel==1) )
            {           
                plot(plotPoints, xlab = xlab, ylab = ylab, log = log, xlim = xLimits, ylim = yLimits, axes = FALSE, 
                     frame.plot = TRUE, col = colourVec[i], pch = pchVec[i], ...) 

                xaxisTicks <- axTicks(1)             
                if (breakCurve) {xaxisTicks <- unique(c(xLimits[1], xaxisTicks))}
                
                yaxisTicks <- axTicks(2)
                axis(2, at=yaxisTicks)
#                xLabels <- as.character(xaxisTicks)
                xLabels <- as.expression(xaxisTicks)
                if (conNameYes) {xLabels[1] <- conName}
                axis(1, at=xaxisTicks, labels=xLabels)
            }
            if (is.null(level) || ((!is.na(matchLevel)) && (matchLevel>1)) )
            {           
                 points(plotPoints, col = colourVec[i], pch = pchVec[i])
            }
#            }
        }
    }
    }
    
    noPlot <- rep(FALSE, numCol)
    for (i in 1:numCol)
    {
#        if (colour) {j <- i} else {j <- 1}  
        if (any(is.na(plotMat[,i]))) {noPlot[i] <- TRUE; next}

        if (is.null(level) || unique(assayNoOld)[i]%in%level)
        {                   
            lines(dosePts, plotMat[,i], lty=ltyVec[i], col=colourVec[i], ...)  # solid only used if there are more curves than types of dashed lines
        }
    }
    
    legendLevels <- as.character(unique(assayNoOld))
    if (!missing(legendText)) 
    {
        lenLT <- length(legendText)
    
        if (lenLT == numAss) {legendLevels <- legendText}
        
        if (lenLT == 1) {legendLevels <- rep(legendText, numAss)}
    }
    levInd <- 1:numAss
    if (!is.null(level)) 
    {
        legendLevels <- level
        levInd <- (1:numAss)[unique(assayNoOld)%in%level]
    } 
    
    
#    ltyIndex <- (1+1:numAss)
    ltyVec[noPlot] <- 0
    
#    if ( (numAss == 1) && (!legend) ) {legend <- FALSE}  # not displaying legend for a single curve unless forced to do so
    if ( legend && !(type == "add") )
    {
        legend(xLimits[2], yLimits[2], legendLevels, lty = ltyVec[levInd], pch = pchVec[levInd], 
               col = colourVec[levInd], bty = "n", xjust = 1, yjust = 1)
    }
    par(las=0)

    retData <- data.frame(dosePts, as.data.frame(plotMat))
#    retData <- as.matrix(retData)
    colnames(retData) <- c(colnames(origData)[1:doseDim], as.character(unique(assayNoOld)))
    
    invisible(retData)
    }    
}
# xyplot(y~x|b, panel=function(x,y,z,subscripts){panel.xyplot(x,y); panel.xyplot(x,z[subscripts], type="l")}, z=z)
