"ED2" <- function(object, percVec, ...)
{

    ## Checking 'percVec' vector ... should be numbers between 0 and 100
    if (any(percVec<=0 | percVec>=100)) {stop("Percentages outside the interval ]0, 100[ not allowed")}
    lenpv <- length(percVec)


    ## Reading data
    curveid <- object$data$curve
    dose <- object$data$dose
    resp <- object$data$resp  # proportions
    total <- object$data$total

    int <- object$data$int
    slope <- object$data$slope
    lower <- object$data$lower
    upper <- object$data$upper

    uniID <- unique(curveid)
    lenuid <- length(uniID)


    lowInd <- !is.null(lower)
    upInd <- !is.null(upper)

    intLevels <- tapply(int, curveid, function(x) x[1])
    sloLevels <- tapply(slope, curveid, function(x) x[1])
    if (lowInd) {lowLevels <- tapply(lower, curveid, function(x) x[1])} else {lowLevels <- rep(NA, lenuid)}
    if (upInd) {upLevels <- tapply(upper, curveid, function(x) x[1])} else {upLevels <- rep(NA, lenuid)}


    ## Obtaining parameter estimates from prediction list
    predList <- object$predict
    parMat <- matrix(0, lenuid, 4)
    parMat[,4] <- rep(1, lenuid)
    for (i in 1:lenuid)
    {
        parMat[i, 1] <- predList$int[intLevels[i]]
        parMat[i, 2] <- predList$slope[sloLevels[i]]
        if (lowInd) {parMat[i, 3] <- predList$low[lowLevels[i]]}
        if (upInd) {parMat[i, 4] <- predList$up[upLevels[i]]}
    }
    parMat[is.na(parMat[, 3]), 3] <- 0
    parMat[is.na(parMat[, 4]), 4] <- 1

#    print(parMat)
#    print(match(parMat[1,],object$fit$par))


    link <- object$link
    logTrans <- object$"logTrans"
    ECp <-
    function(p, par)
    {
#        exp((binomial(link)$linkfun((p-par[3])/(par[4]-par[3])) - par[1])/par[2])
##        retVal <- (binomial(link)$linkfun((p-par[3])/(par[4]-par[3])) - par[1])/par[2]
        retVal <- (binomial(link)$linkfun(p) - par[1])/par[2]
        if (logTrans) {retVal <- exp(retVal)}

        return(retVal)
    }


    if (object$link == "logit")
    {
        derLink <- function(p) {1/(p*(1-p))}
    }
    if (object$link == "probit")
    {
        derLink <- function(p) {1/(dnorm(qnorm(p)))}
    }
    if (object$link == "cloglog")
    {
        derLink <- function(p) {1/((1-p)*log(1-p))}
    }


    ## Obtaining variance-covariance matrix
    varMat <- summary(object)$"varMat"


    ## Calculating effect doses
    dimNames <- rep("", lenuid*lenpv)
    retMat <- matrix(0, lenuid*lenpv, 2)
    for (i in 1:lenuid)
    {
        for (j in 1:lenpv)
        {
            rowIndex <- (i-1)*lenpv+j

            pari <- parMat[i,]
            parIndex <- match(pari,object$fit$par)
            parIndex <- parIndex[!is.na(parIndex)]
#            print(parIndex)

            p <- 1-percVec[j]/100
            ecp <- ECp(p, pari)

#            derVec <- c(-ecp/pari[2], -ecp*log(ecp)/pari[2])
#
#            if (lowInd) {derVec <- c(derVec, (ecp/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[4])/((pari[4]-pari[3])^2)))}
#            if (upInd) {derVec <- c(derVec, -(ecp/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[3])/((pari[4]-pari[3])^2)))}

            if (logTrans)
            {
                derVec <- c(-ecp/pari[2], -ecp*log(ecp)/pari[2])
#                if (lowInd) {derVec <- c(derVec, (ecp/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[4])/((pari[4]-pari[3])^2)))}
                if (lowInd) {derVec <- c(derVec, 0)}
#                if (upInd) {derVec <- c(derVec, -(ecp/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[3])/((pari[4]-pari[3])^2)))}
                if (upInd) {derVec <- c(derVec, 0)}
            } else {
                derVec <- c(-1/pari[2], -ecp/pari[2])
#                if (lowInd) {derVec <- c(derVec, (1/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[4])/((pari[4]-pari[3])^2)))}
                if (lowInd) {derVec <- c(derVec, 0)}
#                if (upInd) {derVec <- c(derVec, -(1/pari[2])*derLink((p-pari[3])/(pari[4]-pari[3]))*((p-pari[3])/((pari[4]-pari[3])^2)))}
                if (upInd) {derVec <- c(derVec, 0)}
            }


            retMat[rowIndex, 1] <- ecp
            retMat[rowIndex, 2] <- sqrt(derVec %*% varMat[parIndex, parIndex] %*% derVec)

            dimNames[rowIndex]<-paste(uniID[i], ":", percVec[j], sep="")
        }

    }
    dimnames(retMat) <- list(dimNames, c("Estimate", "Std. Error"))
    cat("\n")
    printCoefmat(retMat)
    invisible(retMat)
}
