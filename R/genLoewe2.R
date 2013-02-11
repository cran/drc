"genLoewe2" <- function(
fixed = rep(NA, 9), names = c("b1", "b2", "c1", "c2", "d", "e1", "e2", "f1", "f2"), ssfct = NULL)
#fixed = rep(NA, 9), names = c("b1", "b2", "c", "d1", "d2", "e1", "e2", "f1", "f2"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 9
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

        parmMat[, 6] <- parmMat[, 6] / dose[, 1]
        parmMat[, 7] <- parmMat[, 7] / dose[, 2]        

        applyImplicitFct <- function(parmVec)
        {
            if ((!is.finite(parmVec[6])) && (!is.finite(parmVec[7])))
            {
                retVal <- parmVec[5]  # returning the baseline
#                retVal <- parmVec[3]  # returning the baseline
            } else {

                implicitFct <- function(oVal)
                {
                    (1 / parmVec[6]) * ((oVal^(-1 / parmVec[8]) - 1)^(1 / parmVec[1])) + 
                    (1 / parmVec[7]) * ((oVal^(-1 / parmVec[9]) - 1)^(1 / parmVec[2])) - 1
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

                ## Descending formulation
                retVal <- ((parmVec[3] + occuprate * (parmVec[5] - parmVec[3])) / parmVec[6]) * (occuprate^(-1 / parmVec[8]) - 1)^(1 / parmVec[1]) +
                ((parmVec[4] + occuprate * (parmVec[5] - parmVec[4])) / parmVec[7]) * (occuprate^(-1 / parmVec[9]) - 1)^(1 / parmVec[2])
                # parmVec[3] = m1
                # parmVec[4] = m2
                # parmVec[5] = b               
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
            initCD <- findcd(dframe[, 1], dframe[,3])  # keep be feeded into the LL.4 fits!             
            startLL.d1 <- as.vector(coef(drm(dframe[, c(3,1)], fct = LL.4())))
#            print(head(dframe[, c(3,1)]))
#            print(startLL.d1)
            startLL.d2 <- as.vector(coef(drm(dframe[, c(3,2)], fct = LL.4())))
        
            if (startLL.d1[1] < 0)  # condition in terms of standard "drc" parameter "b" (ascending curves)
            {
                initval <- c(startLL.d1, startLL.d2, c(1, 1))[c(1, 5, 3, 7, 2, 4, 8, 9, 10)]
                initval[3:5] <- initCD[c(2, 2, 1)]            
         
            } else {  # positive "b" (descending curves)
                initval <- c(startLL.d1, startLL.d2, c(1, 1))[c(1, 5, 2, 6, 3, 4, 8, 9, 10)] * c(-1, -1, rep(1, 7))
                initval[3:5] <- initCD[c(1, 1, 2)]            
            }       
      
#            initval <- c((llogistic()$ssfct(dframe))[c(1, 1, 2, 2, 3, 4, 4)], 1, 1) * c(-1, -1, rep(1, 7))
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
    name = "genLoewe2",
    text = "Generalized Loewe2", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "genLoewe2"
    invisible(returnList)
}
