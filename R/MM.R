"MM" <- function(
fixed = c(NA, NA, NA), names = c("y0", "Vm", "K"), fctName, fctText)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    ## Fixing parameters
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
#    parmVec1 <- parmVec
#    parmVec2 <- parmVec
    
    ## Defining the non-linear function
    MMfct<-function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        parmMat[, 1] + (parmMat[, 2] - parmMat[, 1]) * dose/(parmMat[, 3] + dose)
    }

    ## Defining self starter function        
    MMssfct<-function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        aPar <- min(y) * 0.95
        bPar <- max(y) * 1.05
    
##        Vm <- max(resp)
#        Vm <- mean(resp[dose == max(dose)])
#        y0 <- min(resp)
##        y0 <- mean(resp[dose == min(dose)])
#        if (length(unique(dose)) == 1) {return((c(NA, NA, y0))[notFixed])}        

        recipX <- 1 / x
        pseudoY <- (bPar - y) / (y - aPar)
        
        cPar <- coef(lm(pseudoY ~ recipX - 1, subset = x > 1e-10))
        return(c(aPar, bPar, cPar)[notFixed])           
        
#        isfi <- is.finite(recipX) & is.finite(recipY)
#        linmodel <- lm(recipY[isfi] ~ recipX[isfi])
##        Vm <- 1/(coef(linmodel)[1])
#        K <- coef(linmodel)[2]*Vm     
    }

    ## Defining names
    MMnames <- names[notFixed]

    ## Defining the ED function
    MMedfct <- function(parm, p, ...)
    {
        parmVec[notFixed] <- parm
    
#        tempVal <- (p/100)*parmVec[1] + (1-p/100)*parmVec[3]
#        EDp <- ((1-p)*(parmVec[1] - parmVec[3])*parmVec[2])/tempVal

        pProp <- p / 100
        tempVal <- (1 - pProp)/(pProp)
        EDp <- parmVec[3] * tempVal  
        EDder <- c(0, 0, tempVal)

#        EDder <- c( ((1-p/100)*parmVec[2]*parmVec[3])/(tempVal^2), 
#                    ((1-p/100)*(parmVec[1] - parmVec[3]))/tempVal,
#                    (-(1-p/100)*parmVec[2]*parmVec[1])/(tempVal^2) )

        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        parmVec[3] * (y - parmVec[1]) / (parmVec[2] - y)
    } 
    
    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        return( cbind(parmMat[, 3], dose, -dose * (parmMat[, 2] - parmMat[ ,1]) / (parmMat[, 3] + dose) ) /
        (parmMat[, 3] + dose) )
    }
    
    ## Defining first derivatives in parameters
    derivx <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
  
        (parmMat[, 2] - parmMat[, 1]) * parmMat[, 3] / ( (parmMat[, 3] + dose)^2 )
    }       

    returnList <- 
    list(fct = MMfct, ssfct = MMssfct, names = MMnames, deriv1 = deriv1, deriv2 = NULL, derivx = derivx,
    edfct = MMedfct, inversion = invfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName), 
    text = ifelse(missing(fctText), "Shifted Michaelis-Menten", fctText), 
    noParm = sum(is.na(fixed)))     
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"MM.2" <-
function(fixed = c(NA, NA), names = c("Vm", "K"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( MM(fixed = c(0, fixed[1:2]), 
    names = c("y0", names[1:2]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Michaelis-Menten")) )
}

"MM.3" <-
function(fixed = c(NA, NA, NA), names = c("y0", "Vm", "K"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( MM(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]])) )
}

#ctl <- data.frame(list(conc=c(0,0,0), rate=c(55,60,65), state=c('ctl','ctl','ctl')))
#pu <- rbind(Puromycin, ctl)
#model1<-multdrc(rate~conc, state, collapse=data.frame(state,state,1), data=pu, fct=MM(fixed=c(NA,NA,NA)))
#summary(model1)
#plot(model1) 
