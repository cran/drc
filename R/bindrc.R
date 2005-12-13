"bindrc" <- function(formula, weights, curve, int, slope, data = NULL, na.action = na.fail, link = "logit", log = TRUE,
                     startVal, lower, upper, fixedLow = NULL, fixedUp = NULL, names = c("a", "b", "c", "d"))
{

    ## Handling the 'formula', 'curveid' and 'data' arguments
    mf <- match.call(expand.dots = FALSE)
    nmf <- names(mf)
    mnmf <- match(c("formula", "weights", "curve", "int", "slope", "data", "na.action", "lower", "upper"), nmf, 0)

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())
    mt <- attr(mf, "terms")

    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
    resp <- model.response(mf, "numeric")
    lenData <- length(resp)

    numTrials <- model.weights(mf)

    curveid2 <- model.extract(mf, "curve")
    intid2 <- model.extract(mf, "int")
    slopeid2 <- model.extract(mf, "slope")    
    lowerid2 <- model.extract(mf, "lower")
    upperid2 <- model.extract(mf, "upper")
    
    
    ## Log-transforming
    if (log) {ldose <- log(dose)} else {ldose <- dose}
        

    ## Finding link
    if (link == "logit" | link == "probit" | link == "cloglog") {famFct <- binomial} else {stop("Link not allowed")}
#    if (link == "probit") {famFct <- binomial(link = "probit")}
#    print(famFct)


    ## Handling 'lowerid', 'upperid', 'fixedLow' and 'fixedUp' arguments

    estIndex <- rep(TRUE, lenData)
    if (!is.null(fixedLow))
    {
        estIndex <- is.na(fixedLow)
    } #else {fixedLow <- rep(0, lenData)}
    if (!is.null(fixedUp))
    {
        estIndex <- is.na(fixedUp)
    } #else {fixedUp <- rep(1, lenData)}

    lowLim <- numeric(0)
    numLL <- 0
    upLim <- numeric(0)
    numUL <- 0
    if (!is.null(lowerid2)) 
    {
        lowLim <- tapply(resp[estIndex], lowerid2[estIndex], function(x) {x[which.min(dose)]})  
        # start values for lower limit(s) ... problem in case of ties
        numLL <- length(lowLim)
    }
    if (!is.null(upperid2)) 
    {
        upLim <- tapply(resp[estIndex], upperid2[estIndex], function(x) {x[which.max(dose)]})  
        # start values for upper limit(s) ... problem in case of ties
        numUL <- length(upLim)
    }    
#    numLL <- length(unique(lowerid))
#    numUL <- length(unique(upperid))
#    lowup <- upper+lower

    lowup <- numLL + numUL
    

    ## Computing start values and constructing design matrix
    if (is.null(curveid2)) 
    {
        curve <- rep(1, lenData)
        intercept <- curve
        slope <- curve
        lenin <- 1
        lensl <- 1
        
        X <- model.matrix(~ldose)
        
        numParm <- lowup + ncol(X)        
        startVal2 <- rep(0, numParm)        

        iv <- is.finite(ldose) & (resp>0) & (resp<1)     
        resp2 <- resp[iv]
        ldose2 <- ldose[iv]     
        startVal2[(lowup+1):numParm] <- coef(lm( famFct(link)$linkfun(resp2) ~ ldose2))
            
    } else {
    
        curve <- curveid2 
        if (!is.null(intid2)) {intercept <- factor(intid2)} else {intercept <- factor(curveid2)}
        if (!is.null(slopeid2)) {slope <- factor(slopeid2)} else {slope <- factor(curveid2)}
        lenin <- length(unique(intercept))
        lensl <- length(unique(slope))
            
        if ((lenin > 1) & (lensl > 1)) 
        {
#            X <- model.matrix(~intercept + slope*ldose)
            X <- cbind(model.matrix(~intercept -1), model.matrix(~slope*ldose - ldose)[,c((lensl+1):(2*lensl))])
        } else {
#            if (lenin > 1) {X <- model.matrix(~intercept - 1 + ldose)}
#            if (lensl > 1) {X <- model.matrix(~slope*ldose - ldose)}
             if (lenin > 1)
             {
                 X <- cbind(model.matrix(~intercept - 1), model.matrix(~ldose - 1))
             }
             if (lensl > 1)
             {
                 X <- model.matrix(~slope*ldose - ldose)[,c(1, (lensl+1):(2*lensl))]
             }
             if ((lenin == 1) & (lensl == 1))
             {
                 X <- model.matrix(~ldose)
             }
        }
#        print(X)
        
        numParm <- lowup + ncol(X)
        startVal2 <- rep(0, numParm)        
        
        iv <- is.finite(ldose) & (resp>0) & (resp<1)
        resp2 <- resp[iv]
        ldose2 <- ldose[iv]
        int2 <- intercept[iv]
        slo2 <- slope[iv]
#        tempSV <- coef(lm( famFct(link)$linkfun(resp2) ~ int2 - 1 + slo2*ldose2 - ldose2))
        tempSV <- coef(lm(resp[iv]~X[iv,] - 1))
#        print(tempSV)
        startVal2[(lowup+1):(numParm)] <- tempSV[!is.na(tempSV)]
    }
    if (lowup > 0) {startVal2[1:lowup] <- c(lowLim, upLim)}
#    print(startVal2)

    X[is.nan(X)] <- 0  # is the always correct? when a contrast is infinite


    ## Assigning lower and upper limits to observation numbers
    asFct <- function(prob, id)
    {
        retVec <- rep(0, lenData)
        
        uniVec <- unique(id)
        for (i in 1:length(uniVec)) {retVec[id%in%c(uniVec[i])] <- prob[i]}
        return(retVec)
    }
#    print(asFct(lowLim, c(lowerid2)))  # coercising into a vector
#    print(asFct(upLim, upperid2))  # coercising into a vector

  
    ## Defining objective function       
    if (!is.null(lowerid2) & !is.null(upperid2))
    {
        optfct <- function(p)
        {           
            p1 <- asFct(p[1:numLL], c(lowerid2)[estIndex])  # coercising into a vector
            p2 <- asFct(p[(numLL + 1):(lowup)], c(upperid2)[estIndex])  # coercising into a vector
            
            prob <- p1 + ( p2 - p1 )*famFct(link)$linkinv(X%*%p[-c(1:lowup)])
    
            return( sum((resp*numTrials)*log(prob/(1-prob))+numTrials*log(1-prob)) )
        }    
    } else if (!is.null(lowerid2))
    {      
        optfct <- function(p)
        {            
#            p1 <- asFct(p[1:numLL], c(lowerid2))  # coercising into vector
            p1 <- asFct(p[1:numLL], c(lowerid2)[estIndex])  # coercising into vector
            
            prob <- p1 + (1 - p1)*famFct(link)$linkinv(X%*%p[-c(1:numLL)])
    
            return( sum((resp*numTrials)*log(prob/(1-prob))+numTrials*log(1-prob)) )
        }

    } else if (!is.null(upperid2))
    {
        p2 <- rep(1, lenData)
        optfct <- function(p)
        {            
#            p2 <- asFct(p[1:numUL], c(upperid2)[estIndex])  # coercising into vector
#            p2 <- rep(1, lenData)
            p2[estIndex] <- asFct(p[1:numUL], c(upperid2))[estIndex]  # coercising into vector                        
                        
            prob <- 0 + (p2 - 0)*famFct(link)$linkinv(X%*%p[-c(1:numUL)])
    
            return( sum((resp*numTrials)*log(prob/(1-prob))+numTrials*log(1-prob)) )
        }    
    } else {
        optfct <- function(p)
        {            
            prob <- 0 + (1 - 0)*famFct(link)$linkinv(X%*%p)
    
            return( sum((resp*numTrials)*log(prob/(1-prob))+numTrials*log(1-prob)) )
        }        
    } 
    if ((!missing(startVal)) && (!(length(startVal2) == length(startVal)))) {stop("Wrong number of start values")}
    if (missing(startVal)) {startVal <- startVal2}
    
#    print(optfct(startVal))  
    
    
    ## Optimising
    fitObject <- try(optim(startVal, optfct, hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)

    if (inherits(fitObject, "try-error")) 
    {
        stop("Convergence failed")
    }


    ## Constructing list for prediction
    if (!is.null(lowerid2)) 
    {
        lowParm <- fitObject$par[1:numLL]
        names(lowParm) <- as.character(unique(lowerid2[estIndex]))
#        print(lowParm)        
    } else {lowParm <- NULL}
    
    if (!is.null(upperid2)) 
    {
        upParm <- fitObject$par[(1+numLL):lowup]
        names(upParm) <- as.character(unique(upperid2[estIndex]))
#        print(upParm)
    } else {upParm <- NULL}

    intParm <- fitObject$par[(lowup+1):(lowup+lenin)]
    names(intParm) <- as.character(unique(intercept))
#    print(intParm)
    
    sloParm <- fitObject$par[(lowup+1+lenin):(lowup+lenin+lensl)]
    names(sloParm) <- as.character(unique(slope))
#    print(sloParm)

    predList <- list(low=lowParm, up=upParm, int=intParm, slope=sloParm)
               
               
    ## Defining labels for the parameters
    parNames <- c(paste(names[1], unique(intercept), sep=":"), paste(names[2], unique(slope), sep=":"))    
    if (!is.null(upperid2)) {parNames <- c(paste(names[4], unique(upperid2[estIndex]), sep=":"), parNames)}
    if (!is.null(lowerid2)) {parNames <- c(paste(names[3], unique(lowerid2[estIndex]), sep=":"), parNames)}
#    print(parNames)
 
 
    ## Fitting saturated model
#    print(X)
    if (is.null(lowerid2) || is.null(upperid2))  # what to do when lower/upper limits are present?
    {
        if (length(unique(curve)) == 1)
        {
           fitSat <- glm(resp~factor(dose), weights=numTrials, family=famFct(link))
        } else {
           fitSat <- glm(resp~factor(curve)*factor(dose), weights=numTrials, family=famFct(link))
        }
    } else {llSat <- NA; dfSat <- NA}
    llSat <- logLik(fitSat)
#    print(fitSat$df.residual)
    dfSat <- fitSat$df.residual    
      
    
    ## Calculating the value of the log-likelihood
    dfFit <- lenData - numParm
    llFit <- fitObject$value + sum(lchoose(numTrials, numTrials*resp))  
    
    
    ## Collecting objects to return
    returnList <- list(fitObject, list(resp=resp, total=numTrials, dose=dose, curve=curve, int=as.character(intercept), slope=as.character(slope), 
                       lower=lowerid2, upper=upperid2, fixedLow=fixedLow, fixedUp=fixedUp), link, X, match.call(), startVal, predList, parNames,
                       c(llSat, dfSat, llFit, dfFit, lenData), log)

    names(returnList) <- c("fit", "data", "link", "design", "call", "startVal", "predict", "parNames", "loglik", "logTrans")
    class(returnList) <- c("drc", "bindrc")
    
    return(returnList)
}
