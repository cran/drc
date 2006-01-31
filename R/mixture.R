"mixture" <- function(formula, curve, collapse, weights, data = NULL, boxcox = FALSE, bcAdd = 0, varPower = FALSE, startVal, fct = l4(), 
                      na.action = na.fail, robust = "mean", type = "continuous", cm = NULL, logDose = NULL, control = mdControl(), 
                      model = "Hewlett", startVal2)
{
    ## Setting na.action option
    options(na.action=deparse(substitute(na.action)))


    ## Checking the 'data' argument 
    if (missing(data)) {stop("'data' argument is required")}


    ## Handling the 'formula', 'curve', 'weights' and 'data' arguments

    mf <- match.call(expand.dots = FALSE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curve", "weights", "data", "na.action"), nmf, 0) 

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())
    mt <- attr(mf, "terms")
        
#    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
#    print(dose)
#    resp <- model.response(mf, "numeric")
#    origResp <- resp  # in case of transformation of the response    
#    lenData <- length(resp)
#    doseDim <- ncol(as.matrix(dose))
#    dimData <- doseDim + 1  # dimension of dose plus 1 dimensional response
    
#    varNames <- names(mf)
#    varNames <- varNames[c(2:dimData,1)]

    ## Handling situation where no 'curve' or 'collapse' argument is given
    assayNo <- model.extract(mf, "curve")
#    collapseId <- model.extract(mf, "collapse")   
    lenData <- length(assayNo)
#    collapseList <- list()


    ## Retrieving weights
    weights <- model.weights(mf)
    if (is.null(weights)) {weights <- rep(1, lenData)}

    
    ## Finding indices for missing values
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}

    ## Removing NAs from the data set
    data2 <- removeMI(data)


    if (is.null(assayNo)) {assayNo <- rep(1, lenData)}

    evalCol <- eval(collapse, envir = data)
    
    if (inherits(evalCol, "formula"))  # (is.list(evalCol)) 
    {
#        if (!(length(evalCol) == 1)) {stop("Two elements should be supplied in 'collapse' list")}
        bCol <- evalCol
#        eCol <- eval(collapse, envir=data)[[2]]
    } else {
        stop("'collapse' argument should be a formula")
    }

    eCol0 <- deparse(substitute(curve))
    eCol <- eval(parse(text=paste("~factor(", eCol0, ")", sep = "")), envir = data)
#    print(eCol)
    
    if (length(fct$names) == 3)  # l3 function 
    {
#        HPfct <- hpfct(fixed = c( NA, 0, NA, NA, NA, NA))
        collapseNew <- list(bCol, ~1, eCol)
        
#        eName <- fct$names[3]
#        noLim <- 1  # number of c and d parameters
    } else if (length(fct$names) == 4)  # l4 function
        {
#            HPfct <- hpfct()
            collapseNew <- list(bCol, ~1, ~1, eCol)
            
#            eName <- fct$names[4]
#            noLim <- 2  # number of c and d parameters  
        } else {stop("Does not work for l5")}    
    assign("collapseNew", collapseNew, envir = .GlobalEnv)


    ## Re-setting na.action
    options(na.action = "na.omit")  # the default


    ## Modifying the control list
    conList <- control
    conList$"noMessage" <- TRUE


    ## Fitting model where ED50 vary freely
    model1 <- multdrc(formula=formula, curve=assayNo, collapse=collapseNew, weights=weights, data=cbind(data2, assayNo, weights), fct=fct, 
                      boxcox=boxcox, bcAdd=bcAdd, varPower=varPower, startVal=startVal, robust=robust, type=type, cm=cm, logDose=logDose, control=conList)
                      # no need for an 'na.action' argument as NAs have been removed above

    rm(collapseNew, envir = .GlobalEnv)


    ## Constructing new collapse argument ... in common for Hewlett and Voelund
    eCol1 <- eval(parse(text=paste("~I(1/(", eCol0, "/100))-1", sep = "")), envir = data)    
    eCol2 <- eval(parse(text=paste("~I(1/(1-", eCol0, "/100))-1", sep = "")), envir = data)
        
#    collapseNew2 <- list(bCol, ~1, ~1, ~I(1/(assayNo/100))-1, ~I(1/(1-assayNo/100))-1, ~1)
#    collapseNew2 <- list(bCol, ~1, ~1, eCol1, eCol2, ~1)

    if (model == "CA")
    {
        if (length(fct$names) == 3)  # l3 function 
        {
            mixtfct <- hewlett(fixed = c( NA, 0, NA, NA, NA, 1))
            collapseNew2 <- list(bCol, ~1, eCol1, eCol2)
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # l4 function
            {
                mixtfct <- hewlett(fixed = c( NA, NA, NA, NA, NA, 1))
                collapseNew2 <- list(bCol, ~1, ~1, eCol1, eCol2)
            
                eName <- fct$names[4]
                noLim <- 2  # number of c and d parameters  
            } else {stop("Does not work for l5")}
        if (missing(startVal2)) {startVal2 <- NULL}
    }
    
    if (model == "Hewlett")
    {
        if (length(fct$names) == 3)  # l3 function 
        {
            mixtfct <- hewlett(fixed = c( NA, 0, NA, NA, NA, NA))
            collapseNew2 <- list(bCol, ~1, eCol1, eCol2, ~1)
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # l4 function
            {
                mixtfct <- hewlett()
                collapseNew2 <- list(bCol, ~1, ~1, eCol1, eCol2, ~1)
            
                eName <- fct$names[4]
                noLim <- 2  # number of c and d parameters  
            } else {stop("Does not work for l5")}
        if (missing(startVal2)) {startVal2 <- 1}
    }

    if (model == "Voelund")
    {
        if (length(fct$names) == 3)  # l3 function 
        {
            mixtfct <- voelund(fixed = c( NA, 0, NA, NA, NA, NA, NA))
            collapseNew2 <- list(bCol, ~1, eCol1, eCol2, ~1, ~1)
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # l4 function
            {
                mixtfct <- voelund()
                collapseNew2 <- list(bCol, ~1, ~1, eCol1, eCol2, ~1, ~1)
            
                eName <- fct$names[4]
                noLim <- 2  # number of c and d parameters  
            } else {stop("Does not work for l5")}
        if (missing(startVal2)) {startVal2 <- c(3, 0.3)}
    }    
        
    assign("collapseNew2", collapseNew2, envir = .GlobalEnv)    


    ## Checking if levels 0 and 100 are present    
    if (all(regexpr("0", as.character(unique(assayNo))) < 0 ))
    {
        stop("Level 0 is missing")   
    }
    if (all(regexpr("100", as.character(unique(assayNo))) < 0 ))
    {
        stop("Level 100 is missing")   
    }
    
    
    ## Constructing starting values 
    sv <- coef(model1)
    
    parNames1 <- model1$"parNames"[[1]]
    parNames2 <- model1$"parNames"[[3]]
    eNames <- as.character(parNames2[regexpr(paste(eName, ":", sep=""), parNames1, fixed = TRUE) > 0])

    pos0 <- match(paste("factor(", eCol0, ")0", sep = ""), eNames)
    pos1 <- match(paste("factor(", eCol0, ")100", sep = ""), eNames)
    if (is.na(pos0)) {pos0 <- 1}  # it is the intercept
    if (is.na(pos1)) {pos1 <- 1}

    noED50 <- length(eNames)
    noB <- length(coef(model1)) - noED50 - noLim

    sv2 <- sv[c(1:(noB + noLim), noB + noLim + pos0, noB + noLim + pos1)]
    sv2[noB+noLim+2] <- sv2[noB+noLim+1] + sv2[noB+noLim+2]
    sv3 <- c(sv2, startVal2)


    ## Fitting Hewlett-Plackett model
    model2 <- multdrc(formula=formula, curve=assayNo, collapse=collapseNew2, weights=weights, data=cbind(data2, assayNo, weights), boxcox=boxcox, bcAdd=bcAdd,
                      varPower=varPower, startVal = sv3, fct = mixtfct, robust=robust, type=type, cm=cm, logDose=logDose, control=conList)
    
    rm(collapseNew2, envir = .GlobalEnv)    

    model1$deviance <- model1$"fit"$"value"
    
    model2$"anova"$"test" <- "F"
    model2$"anova"$"anovaFit" <- model1

    return(model2)
}
