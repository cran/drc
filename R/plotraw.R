"plotraw" <-
function(formula, curveid, data=NULL, na.action=na.fail, trtid, colour=FALSE, conLevel=1e-2, conName="0", trellis=FALSE, log="x", xlab, ylab, ...)
{
    ## Handling formula, curveid, data and na.action arguments
    mf <- match.call(expand.dots = TRUE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curveid", "data", "na.action", "trtid"), nmf, 0) 
#    print(mnmf)

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())
    mt <- attr(mf, "terms")
        
    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
    resp <- model.response(mf, "numeric")
#    origResp <- resp  # in case of transformation of the response    
    lenData <- length(resp)
#    doseDim <- ncol(as.matrix(dose))
#    dimData <- doseDim + 1  # dimension of dose plus 1 dimensional response
    
    
    ## Finding indices for missing values
    missingIndices <- attr(mf,"na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}

    assayNoOld <- model.extract(mf, "curveid")
    if (is.null(assayNoOld)) {assayNoOld <- rep(1, lenData)}
#    print(assayNoOld)
#    print(names(mf))
#    numAss <- length(unique(assayNoOld))


    ## Constructing names for strips
    trt <- model.extract(mf, "trtid")    
    stripNames <- unique(assayNoOld)    
    if (!is.null(trt)) 
#    {
#        stripNames <- unique(assayNoOld)
#    } else {
    {
        stripNames <- paste(stripNames, tapply(as.character(trt), assayNoOld, function(x) {x[1]} ), sep=":")
    }
    if (length(stripNames)==1) {stripNames=""}

    
    ## Extracting names
    varNames <- names(mf)[c(2,1)]  # response first and dose second

    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    }
     
    assayNo <- colConvert(assayNoOld)
    assayLevels <- unique(assayNoOld)
    numAss <- length(assayLevels)
#    numCol <- numAss        
#    varNames <- c("dose", "resp")    
        
    
    ## Handling small dose values
    conNameYes <- FALSE
    if (min(dose) < conLevel)
    {
        smallDoses <- dose<conLevel 
        dose[smallDoses] <- conLevel
        conNameYes <- TRUE
        
#        tickNumber <- trunc(log10(max(dose))-log10(min(dose))+1)
#        print(tickNumber)
#        tickPos <- (10^(as.numeric(formatC(seq(log10(conLevel), log10(max(dose)), length=5), digits=1))))
        
        tickPos <- 10^(seq(round(log10(conLevel))-1, round(log10(max(dose))), by=1))
        
        tickPos[1] <- conLevel
        tickPos <- unique(tickPos)
#        print(tickPos)        
        
        tickLabels <- as.character(tickPos)
        tickLabels[1] <- "0" 
    } else {
        tickPos <- seq(signif(min(dose), 1), signif(max(dose), 1), length=5)
#        print(tickPos)
        tickLabels <- as.character(tickPos)  
#        print(tickLabels)      
    }


    ## Assigning axis and main names
    if (missing(xlab)) {if (varNames[1] == "") {xlab <- "Dose"} else {xlab <- varNames[1]}}
    if (missing(ylab)) {if (varNames[2] == "") {ylab <- "Response"} else {ylab <- varNames[2]}}
 

    ## Plotting trellis
    stripNames <- sort(stripNames)  # to ensure strip names and plotted data coincide
    if (trellis) 
    # & (ncol(data1)==3))
    {
#        print(tickPos)
#        print(tickLabels)
#        print(stripNames)
    
        require(lattice, quietly=TRUE)
        tpfct <- function(MyY, MyDose, MyClass, ...) 
        {
            trellis.device(color=FALSE)
            
            xyplot(MyY~MyDose|as.factor(MyClass), 
                   scales = list(x = list(log = ifelse(log=="x", TRUE, FALSE),
                   at = tickPos,
                   labels = tickLabels)),
                   xlab = xlab,
                   ylab = ylab,
                   panel = function(x, y) {
                   panel.grid()
                   m <- sort.list(x)
                   panel.xyplot(x[m], y[m], type = "p", cex = 1) },
                   strip = strip.custom(factor.levels = stripNames),
                   ...)
        }
        print(tpfct(resp, dose, assayNoOld, ...))
    } else {


    ## Constructing vector of colours
    colourVec <- rep(1, numAss)
    if (is.logical(colour) && colour) {colourVec <- 1:numAss}  
    if (!is.logical(colour) && is.numeric(colour) && length(colour)==numAss) {colourVec <- colour}# else {stop("Argument 'colour' not correct")}


    ## Plotting the observations
    par(las=1)
    for (i in 1:numAss)
    {
#        if (colour) {j <- i} else {j <- 1}
        if (i==1) 
        {
            plot(dose[assayNo==i],resp[assayNo==i], xlab=xlab, ylab=ylab, log=log, xlim=c(min(dose), max(dose)),
                 ylim=c(min(resp), max(resp)), pch=i, axes=FALSE, frame.plot=TRUE, col=colourVec[i], ...) 
            xaxisTicks <- axTicks(1)
            yaxisTicks <- axTicks(2)
            axis(2, at=yaxisTicks)
            xLabels <- as.character(xaxisTicks);
            if (conNameYes) {xLabels[1] <- conName}
            axis(1, at=xaxisTicks, labels=xLabels)
        } else {
            points(dose[assayNo==i], resp[assayNo==i], pch=i, col=colourVec[i])
        }
    }
    legend(max(dose), max(resp), as.character(assayLevels), pch=1:numAss, col=colourVec, bty="n", xjust=1, yjust=1)
    par(las=0)
    }    
}
