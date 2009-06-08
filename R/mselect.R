"mselect" <- function(object, fctList = NULL, nested = FALSE)
{
    if (!is.logical(nested)) {stop("'nested' argument takes only the values: FALSE, TRUE")}
    
    contData <- identical(object$"type", "continuous")
    nestedInd <- 3 + contData + nested
    
#    if (missing(fctList)) 
#    {
#        return(object)
#    } else {
    mc <- match.call()
        
    lenFL <- length(fctList) 
    retMat <- matrix(0, lenFL + 1, 3 + contData + nested)

    retMat[1 ,1] <- logLik(object)
    retMat[1, 2] <- AIC(object)
    retMat[1, 3] <- modelFit(object)[2, 5]
    if (contData) {retMat[1, 4] <- summary(object)$"resVar"}
    if (nested) {retMat[1, nestedInd] <- NA}

#        fctList2 <- list()
#        fctList2[[1]] <- deparse((object$"call"$"fct"))
       
#        retList <- list()
#        retList[[1]] <- object

    fctList2 <- rep("", lenFL + 1)        
#    fctList2[1] <- deparse((object$"call"$"fct"))        
    fctList2[1] <- object$"fct"$"name"
       
    if (!is.null(fctList))
    {
        prevObj <- object    
        for (i in 1:lenFL)
        {
#            tempObj <- update(object, fct = fctList[[i]])  # try(update(object, fct = fctList[[i]]), silent = TRUE)
            tempObj <- try(update(object, fct = fctList[[i]]), silent = TRUE)            

            fctList2[i+1] <- fctList[[i]]$"name"            
            if (!inherits(tempObj, "try-error"))
            {   
#                if (is.null(names(fctList)))
#                {

#                    tempChar <- deparse(mc[[3]][i+1])
#                    fctList2[i+1] <- substr(tempChar, start = 1, stop = nchar(tempChar) - 2)                    

#                    fctList2[i+1] <- fctList[[i]]$"name"

#                } else {
#                    tempChar <- names(fctList)[i]
#                    fctList2[i+1] <- as.character(tempChar)
#                }
                
            
                retMat[i+1, 1] <- logLik(tempObj)
                retMat[i+1, 2] <- AIC(tempObj)
                retMat[i+1, 3] <- modelFit(tempObj)[2, 5]
                
                if (contData)
                {
                    retMat[i+1, 4] <- summary(tempObj)$"resVar"
                }                
                
                if (nested) 
                {
                    retMat[i+1, nestedInd] <- anova(prevObj, tempObj, details = FALSE)[2, 5]  # extracting p-value
                }  
            } else {
                retMat[i+1, ] <- NA
            }
            prevObj <- tempObj
        }
    }
#    }
    rownames(retMat) <- as.vector(unlist(fctList2))
    
    cnames <- c("logLik", "AIC", "Lack of fit")
    if (contData) {cnames <- c(cnames, "Res var")}
    if (nested) {cnames <- c(cnames, "Nested F test")}
    colnames(retMat) <- cnames
    
#    print(retMat)
    return(retMat)
#    return(list(fctList2, retList))
}
