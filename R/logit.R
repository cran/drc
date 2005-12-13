"logit" <-
function(lowerc=c(-Inf, -Inf, -Inf, -Inf), upperc=c(Inf, Inf, Inf, Inf), fixed=c(NA, NA, 0, 1), 
         names=c("a","b","c","d"), useDer=FALSE, logTrans=TRUE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    if (!is.logical(logTrans)) 
    {
        stop("Not logical 'log' argument")
    } else {
        if (logTrans) {logFct <- function(x) {log(x)}} else {logFct <- function(x) {x}}
    }


    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining link function and derivative
    linkfun <- function(p) {log(p/(1-p))}
    
    derLink <- function(p) {1/(p*(1-p))}
    
    
    ## Defining distribution
    dist <- "binomial"
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
#        print(dose)
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
#        print(parmMat)

        ## Establishing continuity at dose=0
        tempVal <- exp(parmMat[,1]+parmMat[,2]*logFct(dose))
        indVec <- !is.finite(tempVal)
        tempVal2 <- tempVal/(1 + tempVal)
        tempVal2[indVec] <- 1
#        print(c(dose[5],parmMat[5,1],parmMat[5,2],logFct(dose)[5]))
 
        parmMat[, 3] + (parmMat[, 4]-parmMat[, 3])*tempVal2
   }


    ## Defining value for control measurements (dose=0)
    confct <- function(drcSign)
    {
        if (drcSign>0) {conPos <- 3} else {conPos <- 4}
        confct2 <- function(parm)
        { 
            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
            parmMat[, notFixed] <- parm
            parmMat[, conPos]
        }
        return(list(pos=conPos, fct=confct2))
    }


    ## Defining flag to indicate if more general ANOVA model
    if (any(is.na(fixed[3:4]))) {anovaYes <- NULL} else {anovaYes <- TRUE}


    ## Defining the self starter function
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0,numParm)

        startVal[4] <- max(resp3)  # upper limit
        startVal[3] <- min(resp3)  # lower limit
               
        if (length(unique(dose2))==1) {return((c(NA, NA, NA, startVal[4]))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]  # improvement using Abbott's correction

        options(warn=-1)
        glmModel <- glm(resp3~logFct(dose3), family=binomial(link="logit"))
        options(warn=0)
        startVal[1:2] <- coef(glmModel)  # the a and b parameter

        return(startVal[notFixed])
    }

   
    ## Defining names
    names <- names[notFixed]


    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL


    ## Limits
    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, p, upper=NULL)  # upper argument not used in 'logit'
    {
        parmVec[notFixed] <- parm

        EDpFct <- function(p, par)
        {
            retVal <- (linkfun(p) - par[1])/par[2]
            if (logTrans) {retVal <- exp(retVal)}

            return(retVal)
        }
        
        tempVal <- 1-p/100
        EDp <- EDpFct(tempVal, parm)

        EDder <- rep(0, numParm)
        if (logTrans)
        {
            EDder[1:2] <- c(-EDp/parm[2], -EDp*log(EDp)/parm[2])
        } else {
            EDder[1:2] <- c(-1/parm[2], -EDp/parm[2])
        }

        return(list(EDp, EDder[notFixed]))
    }


    ## Defining a goodness-of-fit test
    gofTest <- function(object)
    {
        total <- (object$"data")[, 5]
        success <- total*(object$"data")[, 2]
                
        expected <- total*fitted(object) 
        
        ## Pearson's statistic
        c( sum( (success - expected)^2 / (expected*(1 - fitted(object))) ), df.residual(object))    
    }


    ## Defining goodness-of-fit function
    if (all(!is.na(fixed[3:4])))
    {
        anovaTest <- function(formula, ds)
        {
#            count <- resp*weights    
            anovaFit <- glm(formula, family=binomial(link = "logit"), data=ds)
            if (df.residual(anovaFit)>0) 
            {
                return(list(test = "lr", anovaFit = anovaFit))
            } else {
                return(NULL)
            }
        }
    } else {
        anovaTest <- NULL
    }

 
    ## Defining the SI function
    sifct <- NULL


    returnList <- list(dist=dist, fct=fct, ssfct=ssfct, confct=confct, anovaYes=anovaYes, names=names, deriv1=deriv1, deriv2=deriv2, lowerc=lowerLimits, 
                       upperc=upperLimits, edfct=edfct, sifct=sifct)

    class(returnList) <- "logit"
    invisible(returnList)
}
