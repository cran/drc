"MM" <- function(fixed = c(NA, NA, 0), names = c("Vm", "K", "y0"))
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
    
        parmMat[,3] + parmMat[,1] * dose/(parmMat[,2] + dose)
    }
    
    
    ## Defining value for control measurements
    confct <- function(drcSign)
    {
        conPos <- 3
        confct2 <- function(parm)
        { 
            parm[, conPos]
        }
        return(list(pos=conPos, fct=confct2))
    }
    

    ## Defining self starter function        
    MMssfct<-function(data) 
    {       
        dose <- data[, 1]
        resp <- data[, 2]
    
#        Vm <- max(resp)
        Vm <- mean(resp[dose == max(dose)])
        y0 <- min(resp)
#        y0 <- mean(resp[dose == min(dose)])
        if (length(unique(dose)) == 1) {return((c(NA, NA, y0))[notFixed])}        

        recipX <- 1/dose
        recipY <- 1/(resp - y0)
        
        isfi <- is.finite(recipX) & is.finite(recipY)
        linmodel <- lm(recipY[isfi] ~ recipX[isfi])

#        Vm <- 1/(coef(linmodel)[1])
        K <- coef(linmodel)[2]*Vm

        return( c(Vm,K,y0)[notFixed] )        
    }


    ## Defining names
    MMnames <- names[notFixed]


    ## Defining the ED function
    MMedfct <- function(parm, p, ...)
    {
        parmVec[notFixed] <- parm
    
        tempVal <- (p/100)*parmVec[1] + (1-p/100)*parmVec[3]
        EDp <- ((1-p)*(parmVec[1] - parmVec[3])*parmVec[2])/tempVal

        EDder <- c( ((1-p/100)*parmVec[2]*parmVec[3])/(tempVal^2), 
                    ((1-p/100)*(parmVec[1] - parmVec[3]))/tempVal,
                    (-(1-p/100)*parmVec[2]*parmVec[1])/(tempVal^2) )

        return(list(EDp, EDder[notFixed]))
    }
    
#    
#    ## Defining the SI function
#    MMsifct <- function(parm1, parm2, pair)
#    {
#
#        ED1 <- edfct(parm1, pair[1])
#        ED2 <- edfct(parm2, pair[2])
#        SIpair <- ED1[[1]]/ED2[[1]]  # calculating the SI value
#        SIder1 <- ED1[[2]]/ED1[[1]]*SIpair
#        SIder2 <- ED2[[2]]/ED2[[1]]*SIpair
#
#        return(list(SIpair, SIder1, SIder2))
#    }    

    returnList <- 
    list(fct = MMfct, confct = confct, ssfct = MMssfct, names = MMnames, edfct = MMedfct)
    class(returnList) <- "Michaelis-Menten"
    invisible(returnList)
}


#ctl <- data.frame(list(conc=c(0,0,0), rate=c(55,60,65), state=c('ctl','ctl','ctl')))
#pu <- rbind(Puromycin, ctl)
#model1<-multdrc(rate~conc, state, collapse=data.frame(state,state,1), data=pu, fct=MM(fixed=c(NA,NA,NA)))
#summary(model1)
#plot(model1) 
