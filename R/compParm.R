"compParm" <-
function(object, strVal, operator = "/", od = FALSE)
{
    strParm <- object$"parNames"[[1]]
    strParm2 <- object$"parNames"[[3]]
    lenSP <- length(strParm)

#    presentVec <- grep(paste("^+",strVal,":", sep=""), strParm)    
#    presentVec <- grep(paste("^\+", strVal, ":", sep=""), strParm)

    if (inherits(object, "mixdrc")) {sep <- ".{1}"} else {sep <- ":{1}"}
    presentVec <- grep(paste("^", strVal, sep, sep=""), strParm)           
#    presentVec <- grep(paste("^", strVal, ":{1}", sep=""), strParm)    
    lenPV <- length(presentVec)
    if (lenPV<2) {stop("No parameters to compare")}


    ## Extracting information from model fit
    sumObj <- summary(object, od = od)
    parm <- sumObj$"coefficients"
    
    if (inherits(object, "bindrc"))
    {
        varMat <- sumObj$"varMat"
    } else {  # taking different parameterisations into account
#        varMat <- object$"transformation"%*%(sumObj$"varMat")%*%t(object$"transformation")
        varMat <- sumObj$"varMat"
    }
    

    ## Defining comparison function and its derivative
    if (operator=="/")
    {
        hypVal <- 1
        fct <- function(ind) {parm[ind[1]]/parm[ind[2]]}
        dfct <- function(ind){sqrt(c(1/parm[ind[2]],-parm[ind[1]]/(parm[ind[2]]^2))%*%varMat[ind,ind]%*%c(1/parm[ind[2]],-parm[ind[1]]/(parm[ind[2]]^2)))}
#        dfct <- function(ind, vm){sqrt(c(1/parm[ind[2]],-parm[ind[1]]/(parm[ind[2]]^2))%*%vm%*%c(1/parm[ind[2]],-parm[ind[1]]/(parm[ind[2]]^2)))}
    }
    if (operator=="-")
    {
        hypVal <- 0
        fct <- function(ind) {parm[ind[1]]-parm[ind[2]]}
        dfct <- function(ind){sqrt(c(1,-1)%*%varMat[ind,ind]%*%c(1,-1))}
    }


    ## Calculating differences or ratios
    lenRV <- lenPV*(lenPV-1)/2
    cpMat <- matrix(0, lenRV, 4)
    compParm <- rep("", lenRV)
    
    degfree <- df.residual(object)  # $"summary"[6]  # sumObj$loglik[2]
    if (is.null(degfree)) {degfree <- 100}  # ad hoc solution for mixdrc    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (object$"type" == "continuous")
    {
        pFct <- function(x) {pt(x, degfree)}
    } else {
        pFct <- pnorm
    }        
   
    k <- 1
    for (i in 1:lenPV) 
    {
        for (j in 1:lenPV)
        {
            if (j<=i) {next}

#            estRP[k] <- fct(presentVec[c(i,j)])  #parm[i]/parm[j]
            cpMat[k, 1] <- fct(presentVec[c(i,j)])  #parm[i]/parm[j]
            
#            seRP[k] <- dfct(presentVec[c(i,j)])
            cpMat[k, 2] <- dfct(presentVec[c(i,j)])

            tVal <- (cpMat[k, 1] - hypVal)/cpMat[k, 2]
            cpMat[k, 3] <- tVal
            cpMat[k, 4] <- pFct(-abs(tVal)) + (1 - pFct(abs(tVal)))

            compParm[k] <- paste(strParm2[presentVec[c(i,j)]], collapse = operator)
            k <- k+1
        }
    }
    dimnames(cpMat) <- list(compParm, c("Estimate", "Std. Error", "t-value", "p-value"))
    
    cat("\n")
    cat("Comparison of parameter", paste("'",strVal,"'",sep=""), "\n")
    cat("\n")
    printCoefmat(cpMat)
    invisible(cpMat)
}
