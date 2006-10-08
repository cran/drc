"weibull" <-
function(lowerc=c(-Inf, -Inf, -Inf, -Inf), upperc=c(Inf, Inf, Inf, Inf), fixed=c(NA, NA, NA, NA), 
         names=c("b","c","d","e"), scaleDose = TRUE, useDer=FALSE)
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
    
        parmMat[,2] + ( parmMat[,3] - parmMat[,2] ) * exp( - exp( parmMat[,1] *( log(dose) - log(parmMat[,4])) ) )
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
#    anovaYes <- list(bin = !any(is.na(fixed[c(2,3)])) , cont = TRUE)
    binVar <- all(fixed[c(2, 3)]==c(0, 1))
    if (is.na(binVar)) {binVar <- FALSE}
    if (!binVar) {binVar <- NULL}    
    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining the self starter function
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[, 1]
        resp3 <- dataFra[, 2]

        startVal <- rep(0, numParm)

        if (is.na(fixed[2]))
        {
            startVal[2] <- min(resp3)  # the lower bound
        } else {
            startVal[2] <- fixed[2]
        }
        
        if (is.na(fixed[3]))
        {
            startVal[3] <- max(resp3)  # the upper bound
        } else {
            startVal[3] <- fixed[3]
        }
        
        
        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA))[notFixed])}  
        # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

#        loglogTrans <- log(-log((resp3-startVal[2] + 0.001)/(startVal[3]-startVal[2])))  # 0.001 to avoid 0 as argument to log
        
#        loglogTrans <- log(-log(abs(resp3 - startVal[2] - startVal[3]/((pi*pi)^2))/(startVal[3] - startVal[2])))
        loglogTrans <- log(-log((resp3 - startVal[2])/(startVal[3] - startVal[2])))

        isFin <- is.finite(loglogTrans)
        loglogTrans <- loglogTrans[isFin]
        dose3 <- dose3[isFin]

#        print(resp3)
#        print(loglogTrans)
#        print(log(dose3))
        loglogFit <- lm(loglogTrans ~ log(dose3))
        
        if (is.na(fixed[4]))
        {
            startVal[4] <- exp(-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter
        } else {
            startVal[4] <- fixed[4]
        }       
#        startVal[4] <- exp(-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter

        if (is.na(fixed[1]))
        {
            startVal[1] <- coef(loglogFit)[2]  # the b parameter
        } else {
            startVal[1] <- fixed[1]
        }       
#        startVal[1] <- coef(loglogFit)[2]  # the b parameter


        ## Avoiding 0 as start value for lower limit (convergence will fail)
        if ( startVal[2] < 1e-12 ) {startVal[2] <- startVal[3]/10}

#        print(startVal)  
        return(startVal[notFixed])
    }

   
    ## Defining names
    names <- names[notFixed]


    ## Defining parameter to be scaled
    if ( (scaleDose) && (is.na(fixed[4])) ) 
    {
        scaleInd <- sum(is.na(fixed[1:4]))
    } else {
        scaleInd <- NULL
    }
    

    ## Defining derivatives
    ## Defining derivatives
    deriv1 <- function(dose, parm)
              {
                  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
                  parmMat[, notFixed] <- parm

                  t1 <- parmMat[, 3] - parmMat[, 2]
                  t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
                  t3 <- exp(-t2)

                  derMat <- as.matrix(cbind( -t1*t3*t2*(log(dose)-log(parmMat[, 4])), 
                                             1 - t3, 
                                             t3, 
                                             t1*t3*t2*parmMat[, 1]/parmMat[, 4] ))
                  return(derMat[, notFixed])
              }
    deriv2 <- NULL


    ## Limits
    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, p, upper=NULL)  # upper argument not used in 'gompertz'
    {
        parmVec[notFixed] <- parm
    
        tempVal <- log(-log((100-p)/100))
        EDp <- exp(tempVal/parmVec[1] + log(parmVec[4]))

        EDder <- EDp*c(-tempVal/(parmVec[1]^2), 0, 0, 1/parmVec[4])
    
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the SI function
    sifct <- function(parm1, parm2, pair)
    {
        parmVec1[notFixed] <- parm1
        parmVec2[notFixed] <- parm2

        tempVal1 <- log(-log((100-pair[1])/100))
        tempVal2 <- log(-log((100-pair[2])/100))
    
        SIpair <- exp(tempVal1/parmVec1[1] + log(parmVec1[4]))/exp(tempVal2/parmVec2[1] + log(parmVec2[4]))
    
        SIder1 <- SIpair*c(-tempVal1/(parmVec1[1]*parmVec1[1]), 0, 0, 1/parmVec1[4])
        SIder2 <- SIpair*c(tempVal2/(parmVec2[1]*parmVec2[1]), 0, 0, -1/parmVec2[4])
    
        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
    }
    

    returnList <- list(fct=fct, confct=confct, anovaYes=anovaYes, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2,
                       lowerc=lowerLimits, upperc=upperLimits, edfct=edfct, sifct=sifct, scaleInd=scaleInd)

    class(returnList) <- "Weibull"
    invisible(returnList)
}
