"boltzmann" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b","c","d","e","f"))
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    


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

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[3] - resp3)/(resp3 - startVal[2]))  # 0.001 to avoid 0 in the denominator
#        print(logitTrans)

        logitFit <- lm(logitTrans ~ dose3)
        startVal[4] <- -coef(logitFit)[1]/coef(logitFit)[2]  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
               
        return(startVal[notFixed])
    }

   
    ## Defining names
    names <- names[notFixed]


    ## Defininf return list
    returnList <- list(fct=fct, ssfct=ssfct, names=names)

    class(returnList) <- "Boltzmann"
    invisible(returnList)
}
