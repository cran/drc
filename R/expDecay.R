"expDecay" <-
function(fixed = c(NA, NA, NA), names = c("span", "K", "plateau"))
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
    
        parmMat[, 1]*exp(-parmMat[, 2]*x) + parmMat[, 3]
    }
    
#    fct2 <- function(x, parm) {parm[1]*exp(-parm[2]*x) + parm[3]}
#    fct <- fParm(fct2, 2, fixed)
  

    ## Defining self starter function        
    ssfct <- function(data) 
    {       
        x <- data[, 1]
        y <- data[, 2]
    
        if (is.na(fixed[3]))
        {
            plateau <- 0.95*min(y)    
        } else {
            plateau <- fixed[3]
        }
#        ifelse(is.na(fixed[3]), 0.95*min(y))
        
        span <- max(y) - plateau
        
        ## Linear regression in order to determine K        
#        tempY <- log( (y - plateau)/span )
#        K <- -coef(lm(tempY ~ x - 1))

        tempY <- log( (y - plateau) )
        coefVec <- coef(lm(tempY ~ x))
        span <- exp(coefVec[1])
        K <- -coefVec[2]
 
        return( c(span,K,plateau)[notFixed] )        
    }
   

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives
    deriv1 <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        helper1 <- exp(-parmMat[, 2]*x)
        cbind(helper1, -parmMat[, 1]*helper1*x, 1)
    }
    
    deriv2 <- NULL

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, deriv1 = deriv1, deriv2 = deriv2) 
    class(returnList) <- "Exponential decay"
    invisible(returnList)
}
