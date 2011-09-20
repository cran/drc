"genursa" <- function(
fixed = rep(NA, 13), names = c("b1", "b2", "c1", "c2", "d", "e1", "e2", "f1", "f2", "e3", "f3", "b3", "c3"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 13
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    ## Defining bisection algorithm  
    bisec <- function(fu, fuLow, fuHigh)
    {
	      for(k in 1:25) 
        {
	          fuMiddle <- 0.5 * (fuLow + fuHigh)
		        if (fu(fuMiddle) > 0) 
            {
			          fuHigh <- fuMiddle
      		  } else {
			          fuLow <- fuMiddle
      		  }
    	  }    
        list(root = fuMiddle)
    }
    
    ## Defining the nonlinear model function
    fct <- function(dose, parm) 
    {
#        print(parm)
        
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        applyImplicitFct <- function(parmVec)
        {
            if ((!is.finite(parmVec[6])) && (!is.finite(parmVec[7])))
            {
                retVal <- parmVec[5]  # returning the baseline
            } else {

                implicitFct <- function(oVal)
                {
                    (1 / parmVec[6]) * ((oVal^(-1 / parmVec[8]) - 1)^(1 / parmVec[1])) + 
                    (1 / parmVec[7]) * ((oVal^(-1 / parmVec[9]) - 1)^(1 / parmVec[2])) +
                    (parmVec[10] / (parmVec[6] * parmVec[7])) * (oVal^(-1 / parmVec[11]) - 1)^(1 / parmVec[12]) - 1
                }
                bisection <- try(bisec(implicitFct, 0, 1), silent = TRUE)
                if (inherits(bisection, "try-error")) 
                {
                    occuprate <- NA
                } else {
                    occuprate <- bisection$"root"
                }
#                print(c(occuprate, 
#                (1 / parmVec[6]) * (occuprate^(-1 / parmVec[8]) - 1)^(1 / parmVec[1]),
#                (1 / parmVec[7]) * (occuprate^(-1 / parmVec[9]) - 1)^(1 / parmVec[2])
#                ))

#                retVal <- parmVec[5] + occuprate * (((parmVec[3] - parmVec[5]) / parmVec[6]) * (occuprate^(-1 / parmVec[8]) - 1)^(1 / parmVec[1]) +
#                ((parmVec[4] - parmVec[5]) / parmVec[7]) * (occuprate^(-1 / parmVec[9]) - 1)^(1 / parmVec[2]))

                retVal <- ((parmVec[3] + occuprate * (parmVec[5] - parmVec[3])) / parmVec[6]) * (occuprate^(-1 / parmVec[8]) - 1)^(1 / parmVec[1]) +
                ((parmVec[4] + occuprate * (parmVec[5] - parmVec[4])) / parmVec[7]) * (occuprate^(-1 / parmVec[9]) - 1)^(1 / parmVec[2]) +
                (parmVec[13] + occuprate * (parmVec[5] - parmVec[13])) * (parmVec[10] / (parmVec[6] * parmVec[7])) * (occuprate^(-1 / parmVec[11]) - 1)^(1 / parmVec[12])
                
            }       
            return(retVal)
        }      
#        print(apply(parmMat, 1, applyImplicitFct))                  
        apply(parmMat, 1, applyImplicitFct)
    }    

    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            LLstart <- llogistic()$ssfct(dframe)
            initval <- c(LLstart[c(1, 1, 2, 2, 3, 4, 4)], rep(1, 4), LLstart[c(1, 2)]) * c(-1, -1, rep(1, 9), -1, 1)
            return(initval[notFixed])
        }        
    }    
    
    ## Defining names
    names <- names[notFixed]

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL

    ## Returning list with model function information
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct = edfct, sifct = sifct,  # scaleFct = scaleFct,
    name = "ursa2",
    text = "Generalized URSA", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "ursa2"
    invisible(returnList)
}
