"richards" <- function(fixed = c(NA, NA, NA, NA, NA), names = c("c", "d", "delta", "kappa", "gamma"))
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) 
    {
        stop("Not correct 'names' argument")
    }
    if (!(length(fixed) == numParm)) 
    {
        stop("Not correct 'fixed' argument")
    }    


    ## Fixing parameters
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    
    
    ## Defining the non-linear function
    Rfct<-function(dose, parm) 
    {
#        print(parm[1,])
        
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        fac1 <- parmMat[, 2] - parmMat[, 1]
        fac2 <- parmMat[, 3]-1
        fac3 <- (1+fac2*exp(-parmMat[, 4]*(log(dose)-log(parmMat[,5]))))^(1/(-fac2))
        
#        print(parmMat[, 1] + fac1*fac3)
        parmMat[, 1] + fac1*fac3
    }
    

    ## Defining self starter function        
    Rssfct<-function(data)  # following Seber and Wild (1989) pp. 332-337
    {       
        dose <- data[, 1]
        resp <- data[, 2]
    
        startVec <- rep(0, 5)

        startVec[1] <- min(resp)
        startVec[2] <- max(resp)  # mean(resp[dose == max(dose)])

#        midResp <- median(resp)
#        startVec[5] <- (dose[resp == midResp])[1]  # only the first value
        startVec[3] <- 4  # inspired by Seber and Wild (1989) p. 334
        
        yStar <- -log( (((resp-startVec[1])/(startVec[2]-startVec[1]))^(1-startVec[3]) - 1)/(startVec[3]-1) )       
        ysFin <- is.finite(yStar)
#        yStar <- yStar[ysFin]  # removing infinite values
        
        ysModel <- lm(yStar[ysFin] ~ dose[ysFin])
        ysCoef <- coef(ysModel)
        
#        startVec[4] <- ysCoef[2]
        startVec[4] <- -startVec[3]
        startVec[5] <- -ysCoef[1]/ysCoef[2]
        
        return( startVec[notFixed] )        
    }


    ## Defining names
    Rnames <- names[notFixed]

    
    ## Returning functions
    returnList <- list(fct = Rfct, ssfct = Rssfct, names = Rnames,
    name = "richards",
    text = "Richards", 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "Richards"
    invisible(returnList)
}
