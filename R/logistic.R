"logistic" <-
function(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
         names = c("b","c","d","e","f"), scaleDose = TRUE, useDer = FALSE, w = FALSE)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}


    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
#    parmVec1 <- parmVec
#    parmVec2 <- parmVec

           
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        parmMat[,2]+(parmMat[,3]-parmMat[,2])/((1+exp(parmMat[,1]*(log(dose)-log(parmMat[,4]))))^parmMat[,5])
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
#    anovaYes <- list(bin = !any(is.na(fixed[c(2,3,5)])) , cont = TRUE)
    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
    if (is.na(binVar)) {binVar <- FALSE}
    if (!binVar) {binVar <- NULL}
    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining a helper function for the self starter function below 
    "oneTail" <- function(dataFra)
    {
        dose <- dataFra[, 1]
        orderDose <- order(dose)
        dose <- log(dose[orderDose] + 1)
        resp <- dataFra[orderDose, 2]  
        lenData <- length(resp)
        if (lenData < 8) 
        {
            return(NULL)
        }

        ivMin <- dose==min(dose)
        ivMax <- dose==max(dose)    
        slope <- (mean(resp[ivMax]) - mean(resp[ivMin]))/(max(dose) - min(dose))
        inter <- mean(resp[ivMin]) - slope*min(dose)

        uniDose <- unique(dose)
        meanResp <- tapply(resp, dose, mean)

        lineVal <- inter + slope*dose
        diff <- resp - lineVal
        
#        sloVec <- rep(0, lenData - 3)
#        for (i in 1:(lenData-3))
#        {
#            lmModel <- lm(resp[i:(i+3)] ~ dose[i:(i+3)])
#            sloVec[i] <- coef(lmModel)[2]
#        }
#        print(sloVec)
#        stop()

        mp <- NULL
        lmModel <- NULL
        if (sum(diff > 0) < 3) 
        {
            lmModel <- lm(resp[(lenData-3):lenData] ~ dose[(lenData-3):lenData])
            
#            mp[1] <- -sign(coef(lmModel)[2])
#            parme <- -coef(lmModel)[1]/coef(lmModel)[2]
#            if (parme < 0) 
#            {
#                mp[2] <- 1
#            } else {
#                mp[2] <- exp(parme) - 1
#            }
            
#            mp <- c( -sign(coef(lmModel)[2]), abs(-coef(lmModel)[1]/coef(lmModel)[2]) )
#             mp <- mean(exp(dose[(lenData-3):lenData]) - 1)
        }
        if (sum(diff < 0) < 3) 
        {
            lmModel <- lm(resp[1:4] ~ dose[1:4])
#            mp <- c( -sign(coef(lmModel)[2]), abs(-coef(lmModel)[1]/coef(lmModel)[2]) )        
#             mp <- mean(exp(dose[1:4]) - 1)
        }
        if (!is.null(lmModel))
        {
            mp[1] <- -sign(coef(lmModel)[2])
            parme <- -coef(lmModel)[1]/coef(lmModel)[2]
            if (parme < 0) 
            {
                mp[2] <- 1
            } else {
                mp[2] <- exp(parme) - 1
            }
        }
        as.vector(mp)
    }
    
    "weightFct" <- function(lower, upper)
    {
        retFct <- function(x)
        {
            fct1 <- function(x) {1 - (12/(upper-lower)^4)*(x - (upper+lower)/2)^4 }
            fct2 <- function(x) {0.25*exp(-abs(x-lower))}
            fct3 <- function(x) {0.25*exp(-abs(x-upper))}    
    
            indVec <- x>=lower && x<=upper
            indVec1 <- x<lower
            indVec2 <- x>upper
        
            retVec <- rep(0, length(x))
            retVec[indVec] <- fct1(x[indVec])
            if (any(indVec1)) {retVec[indVec1] <- fct2(x[indVec1])}
            if (any(indVec2)) {retVec[indVec2] <- fct3(x[indVec2])}
    
#            if (x>=lower && x<=upper)
#            {
#                1 - (12/(upper-lower)^4)*(x - (upper+lower)/2)^4 
#            } else {
#                0.25*exp(-min(c(abs(x-lower), abs(x-upper))))
#            }
            retVec
        }
        return(retFct)
    }    

    ## Defining the self starter function
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

        startVal[3] <- max(resp3) + 0.001  # the d parameter
#        startVal[3] <- mean(resp3[dose2 == max(dose2)]) + 0.001
        
        startVal[2] <- min(resp3) - 0.001  # the c parameter
#        startVal[2] <- mean(resp3[dose2 == min(dose2)]) + (1e-8)*((max(resp3) - min(resp3))/max(resp3))  

        startVal[5] <- 1  # better choice may be possible!        
#        startVal[!notFixed] <- fixed[!notFixed] 
        
        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[3]-resp3)/(resp3-startVal[2]+0.001))  # 0.001 to avoid 0 in the denominator
        
#        sv <- rep(NA, numParm)
#        sv[!notFixed] <- fixed[!notFixed]
#        print(sv)
#        wf <- weightFct( ifelse(!is.na(sv[2]), sv[2], startVal[2]), startVal[3])
#        print(wf(resp3))        
        if (w) 
        {
            wf <- weightFct(startVal[2], startVal[3])        
            wfVal <- wf(resp3)
        } else {
            wfVal <- rep(1, length(resp3))
        }

#        isfi <- is.finite(dose3)  # removing infinite dose values
        logitFit <- lm(logitTrans ~ log(dose3), weights = wfVal)
        startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
        

        ## Adjusting estimate of ED50 if only one tail is observed        
#        otVal <- oneTail(dataFra)
#        if (!is.null(otVal))
#        {
#            if ( (abs(startVal[4]/otVal[2]) < 1000) && (abs(startVal[4]/otVal[2]) > 0.001) )
#            {
#                startVal[1] <- otVal[1]
#                startVal[4] <- otVal[2]
#            }
#        }

#        retVec <- startVal[notFixed]
#        attr(retVec, "pos") <- (1:numParm)[notFixed]
#        return(retVec)
        
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
    deriv1 <- function(dose, parm)
              {
                  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
                  parmMat[, notFixed] <- parm

                  t1 <- parmMat[, 3] - parmMat[, 2]
                  t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
                  t3 <- (1 + t2)^(2*parmMat[, 5])
                  t4 <- parmMat[, 5]*((1 + t2)^(parmMat[, 5] - 1))
                  t5 <- (1 + t2)^parmMat[, 5]                  

                  derMat <- as.matrix(cbind( -t1*t2*t4*(log(dose)-log(parmMat[, 4]))/t3, 
                                             1 - 1/t5, 
                                             1/t5, 
                                             t1*t2*t4*parmMat[, 1]/parmMat[, 4]/t3, 
                                             -t1*log(1+t2)/t5 ))
                  return(derMat[, notFixed])
              }
    deriv2 <- NULL


    ## Setting the limits
    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, p, upper=NULL)  # upper argument not used in 'logistic'
    {
        parmVec[notFixed] <- parm
    
        tempVal <- log((100-p)/100)
        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])

        EDder <- EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
                       0, 
                       0, 
                       1/parmVec[4], 
                       exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))
    
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the SI function
    sifct <- function(parm1, parm2, pair)
    {
#        parmVec1[notFixed] <- parm1
#        parmVec2[notFixed] <- parm2
#            
#        tempVal1 <- log((100-pair[1])/100)
#        tempVal2 <- log((100-pair[2])/100)
#    
#        SIpair <- (parmVec1[4]/parmVec2[4])*(exp(-tempVal1/parmVec1[5])-1)^(1/parmVec1[1])/((exp(-tempVal2/parmVec2[5])-1)^(1/parmVec2[1]))
#
#        SIder1 <- SIpair*c(-log(exp(-tempVal1/parmVec1[5])-1)/(parmVec1[1]^2),
#                           0,
#                           0,
#                           1/parmVec1[4], 
#                           exp(-tempVal1/parmVec1[5])*tempVal1/(parmVec1[5]^2)*(1/parmVec1[1])*((exp(-tempVal1/parmVec1[5])-1)^(-1))) 
#
#        SIder2 <- SIpair*c(log(exp(-tempVal2/parmVec2[5])-1)/(parmVec2[1]^2),
#                           0,
#                           0,
#                           -1/parmVec2[4],
#                           exp(-tempVal2/parmVec2[5])*tempVal2/(parmVec2[5]^2)*(1/parmVec2[1])*((exp(-tempVal2/parmVec2[5])-1)^(-1))) 
#
#        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))


        ED1 <- edfct(parm1, pair[1])
        ED2 <- edfct(parm2, pair[2])
        SIpair <- ED1[[1]]/ED2[[1]]  # calculating the SI value
        SIder1 <- ED1[[2]]/ED1[[1]]*SIpair
        SIder2 <- ED2[[2]]/ED2[[1]]*SIpair

        return(list(SIpair, SIder1, SIder2))
    }
    
    
    returnList <- list(fct=fct, confct=confct, anovaYes=anovaYes, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2,
                       lowerc=lowerLimits, upperc=upperLimits, edfct=edfct, sifct=sifct, scaleInd=scaleInd)

    class(returnList) <- "logistic"
    invisible(returnList)
}
