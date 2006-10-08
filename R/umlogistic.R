"umlogistic" <- function(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf), 
                         fixed = c(NA, NA, NA, NA, NA), names = c("b","c","d","e","f"), alpha, 
                         scaleDose = TRUE, useDer = FALSE)
{
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}
    
    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec


    ## Defining the function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        fracTemp <-     
        parmMat[, 2] + parmMat[, 3] - ((parmMat[, 3]- parmMat[, 2] + parmMat[, 5]*exp(-1/dose^alpha))/(1+exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))))
    }


    ## Defining self starter function
    ssfct <- function(dataFra)
    {
        firstIV <- logistic()$ssfct(dataFra)
        
        firstIV[1] <- -firstIV[1]
        firstIV[5] <- 0
    
        return(firstIV[notFixed])
    }
    
    
    ## Setting the names of the parameters
    names <- names[notFixed]


    ## Defining parameter to be scaled
    if ( (scaleDose) && (is.na(fixed[4])) ) 
    {
        scaleInd <- sum(is.na(fixed[1:4]))
    } else {
        scaleInd <- NULL
    }


    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL


    ## Setting limits
    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, p, upper = 10000, interval = c(1e-4, 10000))
    {    
        mlogistic(alpha = alpha)$edfct(parm, 100-p, upper, interval) 
    }


    ## Defining the SI function
    sifct <- function(parm1, parm2, pair, upper = 10000, interval = c(1e-4, 10000))
    {
        mlogistic(alpha = alpha)$sifct(parm1, parm2, 100-pair, upper, interval)
    }    


    ## Finding the maximal hormesis
    maxfct <- function(parm, upper)
    {
       retVal <- mlogistic(alpha = alpha)$maxfct(parm, upper)
       retVal[2] <- (parm[2] + parm[3]) - (retVal[2] - parm[2])
       
       return(retVal)
    }


    returnList <- list(fct = fct, confct = confct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
                       lowerc = lowerLimits, upperc = upperLimits, edfct = edfct, sifct = sifct, maxfct = maxfct, 
                       scaleInd = scaleInd, anovaYes = anovaYes)    

    returnList <- list(fct = fct, ssfct = ssfct, names = names, edfct = edfct, sifct = sifct, maxfct = maxfct)
    class(returnList) <- "uml4a"
    invisible(returnList)
}
