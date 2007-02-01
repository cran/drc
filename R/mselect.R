"mselect" <- function(object, fctList = NULL, nested = FALSE)
{
    if (!is.logical(nested)) {stop("'nested' argument takes only the values: FALSE, TRUE")}
#    if (missing(fctList)) 
#    {
#        return(object)
#    } else {
        mc <- match.call()
        
        lenFL <- length(fctList) 
        retMat <- matrix(0, lenFL + 1, 4 + nested)

        retMat[1 ,1] <- logLik(object)
        retMat[1, 2] <- AIC(object)
        retMat[1, 3] <- summary(object)$"resVar"
        retMat[1, 4] <- anova(object)[2, 5]
        if (nested) {retMat[1, 5] <- NA}

#        fctList2 <- list()
#        fctList2[[1]] <- deparse((object$"call"$"fct"))
       
#        retList <- list()
#        retList[[1]] <- object

        fctList2 <- rep("", lenFL + 1)        
        fctList2[1] <- deparse((object$"call"$"fct"))        
    
       
    if (!is.null(fctList))
    {
        prevObj <- object    
        for (i in 1:lenFL)
        {
            tempObj <- update(object, fct=fctList[[i]])  # try(update(object, fct = fctList[[i]]), silent = TRUE)
            
            if (!inherits(tempObj, "try-error"))
            {   
#                if (is.null(names(fctList)))
#                {
                    tempChar <- deparse(mc[[3]][i+1])
                    fctList2[i+1] <- substr(tempChar, start = 1, stop = nchar(tempChar) - 2)
#                } else {
#                    tempChar <- names(fctList)[i]
#                    fctList2[i+1] <- as.character(tempChar)
#                }
                
            
                retMat[i+1, 1] <- logLik(tempObj)
                retMat[i+1, 2] <- AIC(tempObj)
                retMat[i+1, 3] <- summary(tempObj)$"resVar"
                retMat[i+1, 4] <- anova(tempObj)[2, 5]
                
                if (nested) 
                {
                    retMat[i+1, 5] <- anova(prevObj, tempObj, details = FALSE)[2, 5]  # extracting p-value
                }  
            } else {
                retMat[i+1, ] <- NA
            }
            prevObj <- tempObj
        }
    }
#    }
    rownames(retMat) <- as.vector(unlist(fctList2))
    
    cnames <- c("logLik", "AIC", "Res var", "Lack of fit")
    if (nested) {cnames <- c(cnames, "Nested F test")}
    colnames(retMat) <- cnames
    
#    print(retMat)
    return(retMat)
#    return(list(fctList2, retList))
}
