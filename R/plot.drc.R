"plot.drc" <-
function(x, ..., level = NULL, broken = FALSE, col = FALSE, conLevel, conName, 
         grid = 100, legend, legendText, legendPos, type = "average", obs, lty, log = "x", 
         cex, pch, xlab, ylab, xlim, ylim,
         bcontrol = NULL, xt = NULL, xtlab = NULL, yt = NULL, ytlab = NULL, add = FALSE, axes = TRUE)
{    
#    breakCurve <- FALSE
    colour <- col  # replace 'colour' throughout the file

    if (!missing(obs)) {stop("Use 'type' instead of 'obs'")}  # remove in the future

    object <- x
#    zeroEps <- 1e-4 # to avoid log(0) at lower limit of dose range in plot 

#    if (inherits(object, "bindrc"))
#    {
#        plotbin(x, ..., colour = colour, conLevel = conLevel, conName = conName, grid = grid, legend = legend, obs = type, col = col, log = log, 
#                xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
#    } else {


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


    ## Determining logarithmic scales
    logX <- TRUE
    if ( (log == "") || (log == "y") ) 
    {
        logX <- FALSE
    }      

    ## Determining range of dose values
    if (missing(xlim)) 
    {
        xLimits <- c(min(dose), max(dose))
    } else {
        xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
    }

    ## Handling small dose values
    conNameYes <- FALSE
    if ( (xLimits[1] < conLevel) && logX) 
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
        if ( (is.null(logDose)) && (logX) )
        {
           dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), length = grid))
        } else {
           dosePts <- seq(xLimits[1], xLimits[2], length = grid)
        }
    } else {}  # No handling of multi-dimensional dose values


    ## Finding minimum and maximum on response scale
    if (is.null(logDose)) 
    {
        plotMat <- plotFct(dosePts)
    } else {
        plotMat <- plotFct(logDose^(dosePts))
    }
    numCol <- ncol(plotMat)
#    print(dosePts)
#    print(plotMat)

    maxR <- max(resp)
    options(warn = -1)  # suppressing warning in case maximum of NULL is taken
    maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    options(warn=0)  

    if (missing(ylim)) 
    {
        if (missing(xlim))
        {
            yLimits <- c(min(resp), maxPM)
        } else {
            yLimits <- getRange(dose, resp, xLimits)
#            print(yLimits)
        }            
    } else {
        yLimits <- ylim
    }


    ## Cutting away y values (determined by the fitted model) outside the limits
    naPlot <- FALSE  # remove naPlot further down
#    for (i in 1:numCol)
#    {
#        logVec <- !(plotMat[, i] >= yLimits[1]  & plotMat[, i] <= yLimits[2])
#        if ( any(!is.na(logVec)) && any(logVec) )
#        {
#            plotMat[logVec, i] <- NA
#            naPlot <- TRUE
#        }
#    }

    ## Determining presence of legend
    if (missing(legend))
    {
        if (numCol == 1) {legend <- FALSE} else {legend <- TRUE}
    }


    ## Constructing vector of colours
    colourVec <- rep(1, numCol)
    if (is.logical(colour) && colour) 
    {
        colourVec <- 1:numCol
    }  
    if (!is.logical(colour) && length(colour) == numCol) 
    {
        colourVec <- colour
    }
    if (!is.logical(colour) && (!(length(colour) == numCol)) ) 
    {
        colourVec <- rep(colour, numCol)
    }   
     
      
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

    ## Defining cex vector
    if (!missing(cex)) 
    {
        if (length(cex)==1) {cexVec <- rep(cex, numCol)} else {cexVec <- cex}
    } else {
        cexVec <- rep(1, numCol)
    }

    ## Creating a broken axis
    if (broken && logX)
    {
        bList <- list(factor = 2, style = "slash", width = 0.02)
            
        if (!is.null(bcontrol))
        {
            namesBC <- names(bcontrol)
            for (j in 1:length(bcontrol))
            {
                bList[[namesBC[j]]] <- bcontrol[[j]]
            }
        }
        breakStyle <- bList$"style"  # "slash"
        breakWidth <- bList$"width"  # 0.02  # default in axis.break
        clFactor <- bList$"factor"  # 2
                     
        brokenx <- clFactor*conLevel                   
        if ( (log == "x") || (log == "xy") || (log == "yx") )
        {
            xRange <- diff(log(xLimits))
            noplotRange <- exp(c(log(brokenx) - (breakWidth/2)*xRange, log(brokenx) + (breakWidth/2)*xRange))
            ivMid <- dosePts > brokenx            
#            ivLeft <- dosePts < noplotRange[1]
#            ivRight <- dosePts > noplotRange[2]

        } else {
#            xRange <- diff(xLimits)
#            noplotRange <- c(brokenx - (breakWidth/2)*xRange, brokenx + (breakWidth/2)*xRange)
             ivMid <- rep(TRUE, grid)
        }
    } else {
        ivMid <- rep(TRUE, grid)
#        ivLeft <- rep(TRUE, grid)
#        ivRight <- rep(TRUE, grid)
    }


    ## Setting a few graphical parameters
    par(las = 1)
    if (!is.null(logDose)) {if (log == "x") {log <- ""}; if ( (log == "xy") || (log == "yx") ) {log <- "y"}}
    
    
    ## Cutting away original x values outside the limits
    eps1 <- 1e-8
    logVec <- !( (dose < xLimits[1] - eps1) | (dose > xLimits[2] + eps1) )
    dose <- dose[logVec]
    resp <- resp[logVec]
    
    
    ## Plotting    
#    if (!type=="add")
#    if (!add)
#    {
    
    uniAss <- unique(assayNoOld)    
    for (i in 1:numCol)
    {
        plotPoints <- switch(type, "average" = cbind(as.numeric(names(tapVec <- tapply(resp[assayNo == i], 
                                                     dose[assayNo == i], mean))), tapVec),
                                   "none"    = c(max(dosePts) + 10, max(c(maxPM, max(resp))) + 10),
                                   "points"  = cbind(dose[assayNo == i], resp[assayNo == i]),
                                   "all"     = cbind(dose[assayNo == i], resp[assayNo == i]),
                                   "obs"     = cbind(dose[assayNo == i], resp[assayNo == i]))
        
        if (!add)
        {
        if (i==1)
        {            
            if (is.null(level) || uniAss[i]%in%level)
            {           
                plot(plotPoints, xlab = xlab, ylab = ylab, log = log, xlim = xLimits, ylim = yLimits, axes = FALSE, 
                     frame.plot = TRUE, col = colourVec[i], pch = pchVec[i], cex = cexVec[i], ...) 
                             
                yaxisTicks <- axTicks(2)
                yLabels <- TRUE
                if (!is.null(yt)) {yaxisTicks <- yt; yLabels <- yt}
                if (!is.null(ytlab)) {yLabels <- ytlab}                
                if (axes) {axis(2, at = yaxisTicks, labels = yLabels)}
                
                xaxisTicks <- axTicks(1)
#                xLabels <- as.character(xaxisTicks)
                xLabels <- as.expression(xaxisTicks)
#                xLabels <- xaxisTicks
                if (conNameYes) {xLabels[1] <- conName}                                                

                
                ## Avoiding too many tick marks on the x axis
                lenXT <- length(xaxisTicks)
                if (lenXT > 6) 
                {
                    halfLXT <- floor(lenXT/2) - 1
                    chosenInd <- 1 + 2*(0:halfLXT)  # ensuring that control always is present
                    xaxisTicks <- xaxisTicks[chosenInd]
                    xLabels <- xLabels[chosenInd]
                }
                
                if (!is.null(xt)) 
                {
                    if (as.character(xt[1]) == as.character(eval(conName))) 
                    {
                        xaxisTicks <- c(xaxisTicks[1], xt[-1])
#                        xaxisTicks[-1] <- xt[-1]
#                        xLabels[1] <- conName 
#                        xLabels[-1] <- xt[-1]
                        xLabels <- c(conName, xt[-1])                        
                    } else {
                        xaxisTicks <- xt
                        xLabels <- xt
                    }
                }
                if (!is.null(xtlab)) {xLabels <- xtlab}
                if (axes) {axis(1, at = xaxisTicks, labels = xLabels)}

                if (broken) 
                {
                    require(plotrix, quietly = TRUE)
                    axis.break(1, brokenx, style = breakStyle, brw = breakWidth)
                }
            }
        } else {

            matchLevel <- match(unique(assayNoOld)[i], level)
            if ( (!is.null(level)) && (!is.na(matchLevel)) && (matchLevel == 1) )
            {           
                plot(plotPoints, xlab = xlab, ylab = ylab, log = log, xlim = xLimits, ylim = yLimits, axes = FALSE, 
                     frame.plot = TRUE, col = colourVec[i], pch = pchVec[i], cex = cexVec[i], ...) 

                yaxisTicks <- axTicks(2)
                yLabels <- TRUE
                if (!is.null(yt)) {yaxisTicks <- yt; yLabels <- yt}
                if (!is.null(ytlab)) {yLabels <- ytlab}                
                axis(2, at=yaxisTicks)
                
                xaxisTicks <- axTicks(1)                                            
#                xLabels <- as.character(xaxisTicks)
                xLabels <- as.expression(xaxisTicks)
                if (conNameYes) {xLabels[1] <- conName}
                
                if (!is.null(xt)) 
                {
                    if (as.character(xt[1]) == as.character(eval(conName))) 
                    {
                        xaxisTicks <- c(xaxisTicks[1], xt[-1])                    
#                        xaxisTicks[-1] <- xt[-1]
#                        xLabels[1] <- conName 
#                        xLabels[-1] <- xt[-1]
                        xLabels <- c(conName, xt[-1])
                    } else {
                        xaxisTicks <- xt
                        xLabels <- xt
                    }
                }
                if (!is.null(xtlab)) {xLabels <- xtlab}
                axis(1, at=xaxisTicks, labels=xLabels)
            }
            if (is.null(level) || ((!is.na(matchLevel)) && (matchLevel>1)) )
            {           
                 points(plotPoints, col = colourVec[i], pch = pchVec[i])
            }
        }
        
        } else {  # add = TRUE
        
#            if (i == 1) 
#            {
                if (is.null(level) || uniAss[i] %in% level) 
                {
                    points(plotPoints, pch = i,col = colourVec[i], ...)
                }
#            } else {
#                matchLevel <- match(uniAss[i], level)
#                if ((!is.null(level)) && (!is.na(matchLevel)) && (matchLevel == 1)) 
#                {
#                    points(plotPoints, xlab = xlab, ylab = ylab, 
#                           xlim = xLimits, ylim = yLimits, 
#                           pch = i, 
#                           col = colourVec[i], ...)
#                }
#                if (is.null(level) || ((!is.na(matchLevel)) && (matchLevel > 1))) 
#                {
#                    points(plotPoints, pch = i, col = colourVec[i])
#                }
#            }     
        }
    }
#    } else {
#       #added by Xiaoyan--begin 
#       if (type == "add") {
#        for (i in 1:numCol) 
#        {
#            plotPoints <-  cbind(as.numeric(names(tapVec <- tapply(resp[assayNo == i], dose[assayNo == i], mean))), tapVec)
#                
#        }
#        }
#         #added by Xiaoyan--end    
#    }
    
    noPlot <- rep(FALSE, numCol)
    if (!(type == "obs"))
    {
    for (i in 1:numCol)
    {
        if ( any(is.na(plotMat[,i])) & (!naPlot) ) 
        {
            noPlot[i] <- TRUE
            next
        }

        if (is.null(level) || unique(assayNoOld)[i]%in%level)
        {                 
#        print(dosePts[ivMid])
            lines(dosePts[ivMid], plotMat[ivMid, i], lty = ltyVec[i], col = colourVec[i], ...)
#            lines(dosePts[ivLeft], plotMat[ivLeft, i], lty = ltyVec[i], col = colourVec[i])
#            lines(dosePts[ivRight], plotMat[ivRight, i], lty = ltyVec[i], col = colourVec[i])
        }
    }
    }


    ## Defining legend and legend text    
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
    ltyVec[noPlot] <- 0
    
    
    ## Removing line types when lines are not drawn
    if (type == "obs")
    {
        ltyVec[levInd] <- 0
    }
    
    
    ## Removing plot symbol when no points are drawn
    if (type == "none")
    {
        pchVec[levInd] <- NA
    }
    
    
    ## Adding legends
    if (legend)
    {
        ## Defining position of legend
        if (!missing(legendPos))
        {
            if ( (is.numeric(legendPos)) && (length(legendPos) == 2) )
            xlPos <- legendPos[1]
            ylPos <- legendPos[2]
        } else {
            xlPos <- xLimits[2] 
            ylPos <- yLimits[2]
        }
    
        legend(xlPos, ylPos, legendLevels, lty = ltyVec[levInd], pch = pchVec[levInd], 
               col = colourVec[levInd], bty = "n", xjust = 1, yjust = 1)
    }
    par(las=0)

    retData <- data.frame(dosePts, as.data.frame(plotMat))
    colnames(retData) <- c(colnames(origData)[1:doseDim], as.character(unique(assayNoOld)))
    
    invisible(retData)
}


getRange <- function(x, y, xlim)
{
    logVec <- (x>=xlim[1] & x<=xlim[2])
    return(range(y[logVec]))
}
