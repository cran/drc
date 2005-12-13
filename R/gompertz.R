"gompertz" <-
function(lowerc=c(-Inf, -Inf, -Inf, -Inf), upperc=c(Inf, Inf, Inf, Inf), fixed=c(NA, NA, NA, NA), 
         names=c("b","c","d","e"), useDer=FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    
    
    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
    
        parmMat[,2] + ( parmMat[,3] - parmMat[,2] ) * exp( - exp( parmMat[,1] *( log(dose) - parmMat[,4]) ) )
    }


    ## Defining value for control measurements (dose=0)
    confct <- function(drcSign)
    {
        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
        confct2 <- function(parm)
        { 
            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
            parmMat[, notFixed] <- parm
            parmMat[, conPos]
        }
        return(list(pos=conPos, fct=confct2))
    }


    ## Defining flag to indicate if more general ANOVA model
    anovaYes <- TRUE


    ## Defining the self starter function
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)
        startVal[3] <- max(resp3)  # +0.001  # the upper bound
        startVal[2] <- min(resp3)  # -0.001  # the lower bound
#        startVal[!notFixed] <- fixed[!notFixed] 

        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        loglogTrans <- log(-log((resp3-startVal[2] + 0.001)/(startVal[3]-startVal[2])))  # 0.001 to avoid 0 as argument to log
        loglogFit <- lm(loglogTrans~log(dose3))
        startVal[4] <- (-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter
        startVal[1] <- coef(loglogFit)[2]  # the b parameter

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
    edfct <- function(parm, p, upper=NULL)  # upper argument not used in 'gompertz'
    {
        parmVec[notFixed] <- parm
    
        tempVal <- log(-log((100-p)/100))
        EDp <- exp(tempVal/parmVec[1] + parmVec[4])

        EDder <- c(EDp*(-tempVal)/(parmVec[1]^2), 0, 0, EDp)
    
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the SI function
    sifct <- function(parm1, parm2, pair)
    {
        parmVec1[notFixed] <- parm1
        parmVec2[notFixed] <- parm2

        tempVal1 <- log(-log((100-pair[1])/100))
        tempVal2 <- log(-log((100-pair[2])/100))
    
        SIpair <- exp(tempVal1/parmVec1[1] + parmVec1[4])/exp(tempVal2/parmVec2[1] + parmVec2[4])
    
        SIder1 <- SIpair*c(-tempVal1/(parmVec1[1]*parmVec1[1]), 0, 0, 1)
        SIder2 <- SIpair*c(tempVal2/(parmVec2[1]*parmVec2[1]), 0, 0, -1)
    
        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
    }
    

    returnList <- list(fct=fct, confct=confct, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2, lowerc=lowerLimits, upperc=upperLimits, 
                       edfct=edfct, sifct=sifct, anovaYes=anovaYes)

#    returnList <- switch(return, "fct+ss" = list(fct,ssfct,names),
#                                 "fct+ss+der" = list(fct,ssfct,names,deriv1,deriv2),
#                                 "ED" = list(edparm, edfct),
#                                 "SI" = list(siparm, sifct))

    class(returnList) <- "gompertz"
    invisible(returnList)
}
