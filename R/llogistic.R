"llogistic" <- function(
lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf), 
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ss = c("1", "2"))
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
        parm[2] + (parm[3]-parm[2])/((1+(x/parm[4])^parm[1]))^parm[5]
    }

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose)-log(parmMat[, 4]))))^parmMat[, 5])
    }


    ## Defining the self starter function
    if (ss == "2")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

#            startVal <- rep(0, numParm)

#            startVal[3] <- max(resp3) + 0.001  # the d parameter
#            startVal[3] <- ifelse(notFixed[3], 1.05*max(y), fixed[3])
#            startVal[3] <- mean(resp3[dose2 == max(dose2)]) + 0.001
             dVal <- ifelse(notFixed[3], 1.01*max(y), fixed[3])
        
#            startVal[2] <- min(resp3) - 0.001  # the c parameter
#            startVal[2] <- ifelse(notFixed[2], 0.95*min(y), fixed[2])
#            startVal[2] <- mean(resp3[dose2 == min(dose2)]) + (1e-8)*((max(resp3) - min(resp3))/max(resp3))  
            cVal <- ifelse(notFixed[2], 0.99*min(y), fixed[2])

#            if (reps)
#            {
#                cVal0 <- median(y[x == min(x)])
#                dVal0 <- median(y[x == max(x)])            
#                if (cVal0 > dVal0)  # making dVal0 the largest
#                {
#                    tval <- cVal0
#                    cVal0 <- dVal0
#                    dVal0 <- tval
#                }
#        
#                cVal <- ifelse(notFixed[2], 0.95*cVal0, fixed[2])
#                dVal <- ifelse(notFixed[3], 1.05*dVal0, fixed[3])                
#            }  

#            startVal[5] <- 1 
            fVal <- 1  # need not be updated with value in 'fixed[5]'
            # better choice than 1 may be possible!        
        
            if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 


            # Cutting away response values close to d
            indexT1a <- x > 0
            indexT1b <- !(y > 0.95*max(y))
            indexT2 <- c(max((1:length(y))[!(indexT1a & indexT1b)]):length(y))
#            print(indexT2)
#            stop()
        
            x2 <- x[indexT2]
            y2 <- y[indexT2]
#            logitTrans <- log((startVal[3] - y)/(y - startVal[2]))
            logitTrans <- log((dVal - y2)/(y2 - cVal))
            logitFit <- lm(logitTrans ~ log(x2))
#            startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
#            startVal[1] <- coef(logitFit)[2]  # the b parameter
        
            coefVec <- coef(logitFit)
            bVal <- coefVec[2]        
            eVal <- exp(-coefVec[1]/bVal)        
        
        
            return(c(bVal, cVal, dVal, eVal, fVal)[notFixed])
        }
    }
    
    if (ss == "1")
    {
        ssfct <- function(dataFra)
        {
            dose2 <- dataFra[,1]
            resp3 <- dataFra[,2]

            startVal <- rep(0, numParm)

            startVal[3] <- max(resp3) + 0.001  # the d parameter
            startVal[2] <- min(resp3) - 0.001  # the c parameter
            startVal[5] <- 1  # better choice may be possible!        
        
            if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 

            indexT2 <- (dose2>0)
            if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
            dose3 <- dose2[indexT2]
            resp3 <- resp3[indexT2]

            logitTrans <- log((startVal[3]-resp3)/(resp3-startVal[2]+0.001))  
            # 0.001 to avoid 0 in the denominator

            logitFit <- lm(logitTrans ~ log(dose3))
            startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
            startVal[1] <- coef(logitFit)[2]  # the b parameter
        
            return(startVal[notFixed])
        }
    }

   
    ## Defining names
    names <- names[notFixed]
    
    
    ##Defining the first derivatives (in the parameters) 
#    if (useD)
#    {
    ## Constructing a helper function
    xlogx <- function(x, p)
    {
        lv <- (x < 1e-12)
        nlv <- !lv
        
        rv <- rep(0, length(x))
        
        xlv <- x[lv] 
        rv[lv] <- log(xlv^(xlv^p[lv]))
        
        xnlv <- x[nlv]
        rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
    
        rv
    }
        
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
#        t3 <- (1 + t2)^(2*parmMat[, 5])
#        t4 <- parmMat[, 5]*((1 + t2)^(parmMat[, 5] - 1))
        t3 <- parmMat[, 5]*((1 + t2)^(-parmMat[, 5] - 1))
        t5 <- (1 + t2)^parmMat[, 5]                  

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t3,  # *t4/t3, 
               1 - 1/t5, 
               1/t5, 
#               t1*t2*t4*parmMat[, 1]/parmMat[, 4]/t3, 
               t1*t2*t3*parmMat[, 1]/parmMat[, 4], 
               -t1*log(1+t2)/t5 )[, notFixed]
    }
        
    deriv2 <- NULL
#    } else {
#        deriv1 <- NULL
#        deriv2 <- NULL        
#    }


    ##Defining the first derivatives (in the parameters)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        temp1 <- x/parmMat[, 4]          
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4 
    }


    ## Setting the limits
    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

  
    ## The three definitions below are not needed in future ('drm')

    ## Defining parameter to be scaled
    if (is.na(fixed[4]))  #  (scaleDose) && (is.na(fixed[4])) ) 
    {
        scaleInd <- sum(is.na(fixed[1:4]))
    } else {
        scaleInd <- NULL
    }
    ## Defining value for control measurements (dose=0)
    confct <- function(drcSign)
    {
        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
        confct2 <- function(parm)
        { 
            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
            parmMat[, notFixed] <- parm
            parmMat[, conPos]
        }
        return(list(pos=conPos, fct=confct2))
    }
    ## Defining flag to indicate if more general ANOVA model
#    anovaYes <- list(bin = !any(is.na(fixed[c(2,3,5)])) , cont = TRUE)
    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
    if (is.na(binVar)) {binVar <- FALSE}
    if (!binVar) {binVar <- NULL}
    anovaYes <- list(bin = binVar, cont = TRUE)
    


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
    
        tempVal <- log((100-p)/100)
        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])

        EDder <- 
        EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
        0, 0, 1/parmVec[4], 
        exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))

# The next lines are not needed because the lower/upper limits are independent of the parameters
# governing the ED values     
#        if (type == "absolute") 
#        {
#            denom <- (parmVec[3] - parmVec[2])^2
#            EDder <- EDder*c(1, (parmVec[3] - respl)/denom, (respl - parmVec[2])/denom, 1, 1)
#        }
        return(list(EDp, EDder[notFixed]))
    }


#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair)
#    {
#        ED1 <- edfct(parm1, pair[1])
#        ED2 <- edfct(parm2, pair[2])
#        SIpair <- ED1[[1]]/ED2[[1]]  # calculating the SI value
#        SIder1 <- ED1[[2]]/ED1[[1]]*SIpair
#        SIder2 <- ED2[[2]]/ED2[[1]]*SIpair
#
#        return(list(SIpair, SIder1, SIder2))
#    }
    
    
    ## Identifying parameters that are on the same scale as x and y (not used)
    if (is.na(fixed[4]))
    {
        sxInd <- sum(is.na(fixed[1:4]))  # sxInd <- c(4)
    } else {
        sxInd <- NULL
    }
    if ( (is.na(fixed[2])) || (is.na(fixed[3])) )
    {
        syInd <- c(sum(is.na(fixed[1:2])), sum(is.na(fixed[1:3])))  # syInd <- c(2, 3)
        if (syInd[1] == 0) {syInd <- syInd[2]}
        if (syInd[2] == 0) {syInd <- syInd[1]}
    } else {
        syInd <- NULL
    }
    
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, bfct = bfct, sxInd = sxInd, syInd = syInd,
    scaleInd = scaleInd, confct=confct, anovaYes=anovaYes, lowerc=lowerLimits, upperc=upperLimits)
    # the last line is not needed in the future ('drm')
    class(returnList) <- "llogistic"
    invisible(returnList)
}

"LL.2" <-
function(fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, 1, fixed[2], 1), names = c(names[1], "c", "d", names[2], "f"), ...) )
}

l2 <- LL.2

"LL.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f"), ...) )
}

l3 <- LL.3

"LL.3u" <-
function(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1:2], 1, fixed[3], 1), names = c(names[1:2], "d", names[3], "f"), ...) )
}

l3u <- LL.3u

"LL.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)  #, reps = FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed, 1), names = c(names, "f"), ...) )  # , reps = reps) )
}

l4 <- LL.4

"LL.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return(llogistic(fixed = fixed, names = names, ...))
}

l5 <- LL.5

