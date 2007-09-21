"llogistic2" <- function(
lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf), 
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ss = c("1", "2", "3"), ssfct = NULL)
{
    ## Matching 'adjust' argument
    ss <- match.arg(ss)
    
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the basic non-linear function
    bfct <- function(x, parm)
    {
        parm[2] + (parm[3]-parm[2])/((1+(x/exp(parm[4]))^parm[1]))^parm[5]
    }

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose)-parmMat[, 4])))^parmMat[, 5])
    }


    ## Defining the self starter function
    
    ## Defining self starter function based on supplied function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        
    ## Version 1 (default)    
    if (ss == "1")
    {
        ssfct <- function(dataFra)
        {
            x <- dataFra[, 1]
            y <- dataFra[, 2]

            startVal <- rep(0, numParm)

            startVal[3] <- max(y) + 0.001  # the d parameter
            startVal[2] <- min(y) - 0.001  # the c parameter
            startVal[5] <- 1  # better choice may be possible!        
        
            if (length(unique(x))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 

            indexT2 <- (x > 0)
#            if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
            x2 <- x[indexT2]
            y2 <- y[indexT2]

            startVal[c(1,4)] <- find.be2(x2, y2, startVal[2] - 0.001, startVal[3])
            startVal[4] <- log(startVal[4])
            # 0.001 to avoid 0 in the denominator

            return(startVal[notFixed])
        }
    }

    ## Version 2
    if (ss == "2")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

            dVal <- ifelse(notFixed[3], 1.01*max(y), fixed[3])        
            cVal <- ifelse(notFixed[2], 0.99*min(y), fixed[2])

            fVal <- 1  # need not be updated with value in 'fixed[5]'
            # better choice than 1 may be possible!        
        
            if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 

            # Cutting away response values close to d
            indexT1a <- x > 0
            x2 <- x[indexT1a]
            y2 <- y[indexT1a]
            
            beVec <- find.be2(x2, y2, cVal, dVal)
            bVal <- beVec[1]
            eVal <- log(beVec[2])
            
            return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        }
    }

    ## Version 3
    if (ss == "3")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

            cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2])
            dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3])
            fVal <- 1  # need not be updated with value in 'fixed[5]'
        
            if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 
           
            beVec <- find.be1(x, y, cVal, dVal)
            bVal <- beVec[1]
            eVal <- log(beVec[2])
            
            return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        }
    }
    }

    ## Finding b and e based on stepwise increments
    find.be1 <- function(x, y, c, d)
    {
        unix <- unique(x)
        uniy <- tapply(y, x, mean)
        lenx <- length(unix)
        
        j <- 2
        for (i in 2:lenx)
        {
            crit1 <- (uniy[i] > (d + c)/2) && (uniy[i-1] < (d + c)/2)
            crit2 <- (uniy[i] < (d + c)/2) && (uniy[i-1] > (d + c)/2)
            if (crit1 || crit2) break
            j <- j + 1
        }
        eVal <- (unix[j] + unix[j-1])/2
        bVal <- -sign(uniy[j] - uniy[j-1])  # -(uniy[j] - uniy[j-1]) / (unix[j] - unix[j-1])
        return(as.vector(c(bVal, eVal)))  
    }
    
    ## Finding b and e based on linear regression
    find.be2 <- function(x, y, c, d)
    {
        logitTrans <- log((d - y)/(y - c))  

        lmFit <- lm(logitTrans ~ log(x))

        coefVec <- coef(lmFit)
        bVal <- coefVec[2]        
        eVal <- exp(-coefVec[1]/bVal)    

        return(as.vector(c(bVal, eVal)))
    }

   
    ## Defining names
    names <- names[notFixed]
    
    
    ##Defining the first derivatives (in the parameters)        
    deriv1 <- NULL 
    deriv2 <- NULL


    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        temp1 <- exp(parmMat[, 1]*(log(x) - parmMat[, 4]))  # x/parmMat[, 4]          
        temp2 <- 1 + temp1
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*temp1*parmMat[, 1]/x
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4 
    }


    ## Setting the limits
    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
        if (type == "absolute") 
        {
            p <- 100*((parmVec[3] - respl)/(parmVec[3] - parmVec[2]))
        } else {  
            p <- respl
        }
        if ( (parmVec[1] < 0) && (reference == "control") )
        {
            p <- 100 - p
        }
    
        tempVal1 <- 100/(100-p)
        tempVal2 <- log(tempVal1^(1/parmVec[5]) - 1)
        lEDp <- parmVec[4] + tempVal2 / parmVec[1]

        lEDder <- 
        c(-tempVal2/(parmVec[1]^2), 
        0, 0, 1, 
        tempVal1^(1/parmVec[5]-1)/(parmVec[1]*parmVec[5]*(tempVal1^(1/parmVec[5]-1))))

        return(list(lEDp, lEDder[notFixed]))
    }
 
    
    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))
    }    
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, bfct = bfct, inversion = invfct)
    class(returnList) <- "llogistic"
    invisible(returnList)
}

"LL2.2" <-
function(fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1], 0, 1, fixed[2], 1), names = c(names[1], "c", "d", names[2], "f"), ...) )
}

l2 <- LL.2

"LL2.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f"), ...) )
}

l3 <- LL.3

"LL2.3u" <-
function(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1:2], 1, fixed[3], 1), names = c(names[1:2], "d", names[3], "f"), ...) )
}

l3u <- LL.3u

"LL2.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)  #, reps = FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed, 1), names = c(names, "f"), ...) )  # , reps = reps) )
}

l4 <- LL.4

"LL2.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return(llogistic2(fixed = fixed, names = names, ...))
}

l5 <- LL.5

