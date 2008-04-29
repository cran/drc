"expDecay" <-
function(fixed = c(NA, NA, NA), names = c("init", "plateau", "k"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    
    ## Defining the non-linear function
    fct <- function(x, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        parmMat[, 1] + (parmMat[, 2] - parmMat[, 1]) * exp(-parmMat[, 3] * x) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        init <- min(y) * 0.95
        plateau <- max(y) * 1.05
        
        ## Linear regression on pseudo y values through origin
        ##  to determine the parameter c
        pseudoY <- log( (y - init) / (plateau - init) )
        kpar <- coef(lm(pseudoY ~ I(-x) - 1))
        return(c(init, plateau, kpar)[notFixed])

#        if (is.na(fixed[3]))
#        {
#            plateau <- 0.95*min(y)    
#        } else {
#            plateau <- fixed[3]
#        }
##        ifelse(is.na(fixed[3]), 0.95*min(y))
#        
#        span <- max(y) - plateau
#        
#        ## Linear regression in order to determine K        
##        tempY <- log( (y - plateau)/span )
##        K <- -coef(lm(tempY ~ x - 1))
#
#        tempY <- log( (y - plateau) )
#        coefVec <- coef(lm(tempY ~ x))
#        span <- exp(coefVec[1])
#        K <- -coefVec[2]
# 
#        return( c(span,K,plateau)[notFixed] )        
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives
    deriv1 <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        tempVal <- exp(-parmMat[, 3] * x)
        cbind(1 - tempVal, tempVal, (parmMat[, 1] - parmMat[, 2]) * tempVal * x)
    }
    
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        ## Creating the full parameter vector
        ##  containing both estimated and fixed parameters
        parmVec[notFixed] <- parm
        
        ## Converting to relative scale
        if (type == "absolute") 
        {
            p <- 100*((parmVec[2] - respl)/(parmVec[2] - parmVec[1]))
        } else {  
            p <- respl
        }
        
        ## Calculating ED value relative to the control
        if (reference == "control")
        {
            p <- 100 - p
        }
        pProp <- p / 100    
        EDp <- -log(1 - pProp)/parmVec[3]

        EDder <- 
        c(0, 0, log(1 - pProp) / (parmVec[3]^2)) 

        return(list(EDp, EDder[notFixed]))
    }
       
    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        -log( (y - parmVec[1]) / (parmVec[2] - parmVec[1]) ) / parmVec[3]
    }    
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, deriv1 = deriv1, deriv2 = deriv2,
    edfct = edfct, inversion = invfct, 
    name = "expDecay",
    text = "Exponential decay",    
    noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}
