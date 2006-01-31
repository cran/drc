"summary.mixdrc" <-
function(object, ...)
{
    objectCopy <- object
    class(objectCopy) <- c("nlme", "lme")
        
    sumObj <- summary(objectCopy)
            
    resultMat <- as.matrix(sumObj$tTable[, c(1,2,4,5)]) 

    varComp <- matrix(as.numeric(VarCorr(objectCopy)[,1]))
    colnames(varComp) <- "Variance"
    rownames(varComp) <- rownames(VarCorr(objectCopy))

    ll <- logLik(objectCopy)
    loglik <- ll[1] 
    degfre <- sumObj$dims$N - attr(ll, "df")

    estimates <-  as.vector(sumObj$coefficients$fixed)
    parNames <- rownames(resultMat)
    varMat <- sumObj$varFix
        
    ## Defining return list
    retList <- list(varComp, varMat, resultMat, c(loglik, degfre), parNames, "mixed logistic")

    names(retList) <- c("varComp", "varMat", "estimates", "loglik", "parNames", "class") 
    class(retList) <- c("summary.drc")
    return(retList)
}
