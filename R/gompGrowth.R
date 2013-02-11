"gompGrowth.1" <-
function(fixed = c(NA, NA, NA), names = c("c", "m", "plateau"))
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
    
        c <- parmMat[, 1]; m <- parmMat[, 2]; a <- parmMat[, 3] 
        a * exp( - (m/c) * exp (-c * x)) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        plateau <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( - log( y / plateau ) )
        coefs <- coef( lm(pseudoY ~ x) )
        k <- coefs[1]; c <- - coefs[2]
        b <- exp(k) 
        m <- b * c
	      return(c(c, m, plateau)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Gompertz Growth Model"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"gompGrowth.2" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "plateau"))
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
    
        a <- parmMat[, 3]; b <- parmMat[, 1]; c <- parmMat[, 2] 
        a * exp( - exp (b * (c - x)))  
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        a <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( log( a / (y +0.0001) ) )
        coefs <- coef( lm(pseudoY ~ x) )

        k <- coefs[1]
        b <- - coefs[2]
        c <- k/b
        
        return(c(b, c, a)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Gompertz Growth Model 2"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

 "gompGrowth.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "c", "plateau"))
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
    
        a <- parmMat[, 3]; b <- parmMat[, 1]; c <- parmMat[, 2] 
        a * exp( - b * exp (-c * x))  
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        a <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( - log( y / a ) )
        coefs <- coef( lm(pseudoY ~ x) )

        k <- coefs[1]
        c <- - coefs[2]
        b <- exp(k)
        
        return(c(b, c, a)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Gompertz Growth Model 3"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}
