"predict.drc" <- function(object, newdata, se.fit = TRUE, interval = c("none", "confidence", "prediction"), 
level = 0.95, na.action = na.pass, od = FALSE, ...)
{
    ## Checking arguments
    interval <- match.arg(interval)

    if (missing(newdata)) 
    {
        predValues <- fitted(object)
        newdata <- data.frame(object$data[, 1], object$data[, 4])
    } 
#    else {
    if (ncol(newdata) > 2) {stop("More than 2 variables in 'newdata' argument")}
    
    
    ## Defining dose values    
    doseVec <- newdata[ ,1]


    ## Constructing matrix of parameter estimates
    parmMat <- object$"parmMat"        
    parmNames <- colnames(parmMat)
    lenCN <- length(parmNames)
    indVec <- 1:lenCN
    names(indVec) <- parmNames
       
    if (lenCN > 1)
    {
        indVec <- indVec[as.character(newdata[, 2])]
        pm <- parmMat[, as.character(newdata[, 2])]  # 'as.character' used to suppress factor levels            
    } else {
        lenDV <- length(doseVec)
        indVec <- rep(1, lenDV)
        pm <- matrix(parmMat[, 1], length(parmMat[, 1]), lenDV)
    }


    ## Checking for NAs in matrix of parameter estimates
#    ncPM <- ncol(parmMat) 
#    naVec <- rep(FALSE, ncPM)
    naVec <- rep(FALSE, lenCN)
    for (i in 1:lenCN)
    {
        naVec[i] <- any(is.na(parmMat[, i]))
    }
    parmMat <- parmMat[, !naVec, drop=FALSE] 


    ## Constructing variance-covariance matrix
    sumObj <- summary(object, od = od)
    varMat <- object$"transformation" %*% sumObj$"varMat" %*% t(object$"transformation") 
    ncPM2 <- ncol(parmMat)
    nrPM2 <- nrow(parmMat)
    indexMat <- matrix(1:(nrPM2*ncPM2), nrPM2, ncPM2, byrow = TRUE)   
    
    
    ## Calculating predicted values  
    indexVec <- as.vector(indVec)    
    lenIV <- length(indexVec)    
    retMat <- matrix(0, lenIV, 4)
    colnames(retMat) <- c("Fit", "SE", "Lower", "Upper")
    retMat[, 1] <- object$"fct"$"fct"(doseVec, t(pm))


    ## Calculating the quantile to be used in the confidence intervals
    if (!(interval == "none"))
    {    
        tquan <- qt(1 - (1 - level)/2, df.residual(object))   
    }
    
    
    ## Calculating standard errors and/or confidence intervals
    if (se.fit || (!interval == "none") )
    {
        rowIndex <- 1    
        for (i in indexVec)
        {
            parmInd <- indexMat[, i]
            varCov <- varMat[parmInd, parmInd]
            parmChosen <- t(parmMat[, i, drop = FALSE])

#            print(parmChosen)
            dfEval <- object$"fct"$"deriv1"(doseVec[rowIndex], parmChosen)
#            print(dfEval)
            varVal <- dfEval%*%varCov%*%dfEval
#            print(varVal)
            retMat[rowIndex, 2] <- sqrt(varVal)  
            
            if (interval == "confidence")
            {
                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal)
                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal)            
            }
            if (interval == "prediction")
            {
                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObj$"resVar")
                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObj$"resVar")                        
            }          
            rowIndex <- rowIndex + 1        
        }
    }
    keepInd <- 1
    if (se.fit) {keepInd <- c( keepInd, 2)}
    if (!interval == "none") {keepInd <- c( keepInd, 3, 4)}
    
    return(retMat[, keepInd])
}


