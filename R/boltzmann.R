"boltzmann" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    


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
        
        parmMat[,2]+(parmMat[,3]-parmMat[,2])/((1+exp(parmMat[,1]*(dose-parmMat[,4])))^parmMat[,5])
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

        logitTrans <- log((startVal[3] - resp3)/(resp3 - startVal[2]))  # 0.001 to avoid 0 in the denominator
#        print(logitTrans)

        logitFit <- lm(logitTrans ~ dose2)
        startVal[4] <- -coef(logitFit)[1]/coef(logitFit)[2]  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
               
        return(startVal[notFixed])
    }

   
    ## Defining names
    names <- names[notFixed]


    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(dose - parmMat[, 4]))
        t3 <- (1 + t2)^(2*parmMat[, 5])
        t4 <- parmMat[, 5]*((1 + t2)^(-parmMat[, 5] - 1))
        t5 <- (1 + t2)^(-parmMat[, 5])                  

        cbind( -t1*t2*t4*(dose-parmMat[ ,1]), 
               1 - 1/t5, 
               1/t5, 
               t1*t2*t4*parmMat[, 1], 
               -t1*log(1+t2)*t5 )[, notFixed]
    }
        
    deriv2 <- NULL
    
    
    ##Defining the first derivatives (in x)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        temp1 <- exp(parmMat[, 1]*(x - parmMat[, 4]))       
        
        (-parmMat[, 5]*(parmMat[, 3] - parmMat[, 2])*temp1*parmMat[, 1])/((1 + temp1)^(parmMat[, 5] + 1))
    }


    ## Defining the ED function
    edfct <- function(parm, p, ...)
    {
        parmVec[notFixed] <- parm
    
        tempVal <- -log((100-p)/100)
        EDp <- parmVec[4] + log(exp(tempVal/parmVec[5])-1)/parmVec[1]

        EDder <- c(-log(exp(tempVal/parmVec[5])-1)/(parmVec[1]^2), 
                   0, 
                   0, 
                   1, 
                   -exp(tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(tempVal/parmVec[5])-1)^(-1)))
    
        return(list(EDp, EDder[notFixed]))
    }


#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair)
#    {
#        ED1 <- edfct(parm1, pair[1])
#        ED2 <- edfct(parm2, pair[2])
#        SIpair <- ED1[[1]] - ED2[[1]]  # SI value on log scale
#        SIder1 <- ED1[[2]]
#        SIder2 <- ED2[[2]]
#
##        SIpair <- ED1[[1]]/ED2[[1]]  # calculating the SI value
##        SIder1 <- ED1[[2]]/ED1[[1]]*SIpair
##        SIder2 <- ED2[[2]]/ED2[[1]]*SIpair
#
#        return(list(SIpair, SIder1, SIder2))
#    }

    ## Defining return list
    returnList <- list(fct = fct, ssfct = ssfct, names = names, 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct)

    class(returnList) <- "Boltzmann"
    invisible(returnList)
}


"B.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( boltzmann(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f")) )
}

b3 <- B.3

"B.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( boltzmann(fixed = c(fixed, 1), names = c(names, "f")) )
}

b4 <- B.4

"B.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))
{
    return( boltzmann(fixed = fixed, names = names) )
}

b5 <- B.5
