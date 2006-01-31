"sumbin" <- function(object, ...)
{
    ## Producing a summary of model fit

    fit <- object$fit
    numParm <- length(fit$par)
    parNames <- object$parNames

    ## Calculating variance-covariance matrix from Hessian
    varMat <- solve(-fit$hessian)

    resultMat <- matrix(0, numParm, 4, dimnames=list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))
    resultMat[,1] <- fit$par

    resultMat[,2] <- sqrt(diag(varMat))

    resultMat[,3] <- resultMat[,1]/resultMat[,2]
    resultMat[,4] <- (pnorm(-abs(resultMat[,3])))+(1-pnorm(abs(resultMat[,3])))


    ## Calculating the value of the log-likelihood
    degfre <- length(object$data$resp) - numParm
    loglik <- fit$value + sum(lchoose(object$data$total, object$data$total*object$data$resp))

    classObj <- class(object)
    if (inherits(object, "list")) {modelClass <- classObj[length(classObj)-1]} else {modelClass <- classObj[length(classObj)]}

    sumObj <- list(varMat, resultMat, c(loglik, degfre), modelClass)
    names(sumObj) <- c("varMat", "estimates", "loglik", "class")
    class(sumObj) <- c("summary.drc")
    return(sumObj)
}
