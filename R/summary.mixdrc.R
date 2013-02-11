"summary.mixdrc" <-
function(object, ...)
{
    objectCopy <- object
    class(objectCopy) <- c("nlme", "lme")
        
    sumObj <- summary(objectCopy)
            
    resultMat <- as.matrix(sumObj$tTable[, c(1,2,4,5)]) 

    varComp <- as.numeric(VarCorr(objectCopy)[,1])

    varComp <- matrix(varComp[!is.na(varComp)])
    colnames(varComp) <- "Variance"
    
    rn <- rownames(VarCorr(objectCopy))   
    rn <- rn[!(rn == "")]
    lenRN <- length(rn)
    if (lenRN == 1) 
    {
        titAt <- attr(VarCorr(objectCopy), "title")
#        esPos <- regexpr("=", titAt)
#        rn <- c( substr(titAt, 1, esPos-1), rn)
        rn <- c( titAt, rn )
    }
    
    for (i in 1:lenRN)
    {
        esPos <- regexpr("=", rn[i])
        if (esPos > 0) {rn[i] <- substr(rn[i], 1, esPos-1)}    
    }
    
    rownames(varComp) <- rn

    ll <- logLik(objectCopy)
    loglik <- ll[1] 
    degfre <- sumObj$dims$N - attr(ll, "df")

#    estimates <-  as.vector(sumObj$coefficients$fixed)
    parNames <- rownames(resultMat)
    varMat <- sumObj$varFix       
    fctName <- paste("mixed", deparse(object$base$call$fct))
        
    ## Defining return list
    retList <- list(varComp, varMat, resultMat, c(loglik, degfre), parNames, "mixed logistic", fctName)

    names(retList) <- c("varComp", "varMat", "coefficients", "loglik", "parNames", "class", "fctName") 
    class(retList) <- c("summary.drc")
    return(retList)
}
