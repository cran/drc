"maED" <- function(object, fctList = NULL, respLev, interval = c("none", "buckland", "kang"), 
level = 0.95, display = TRUE, na.rm = FALSE, extended = FALSE)
{
    interval <- match.arg(interval)

    msMat <- do.call("mselect", list(object = object, fctList = fctList))
       
    expVec <- as.vector(exp(-msMat[, 2] / 2))
    wVec <- expVec / sum(expVec, na.rm = na.rm)  
    # maybe better "combined" na.rm approach for edEst and wVec
    
#    ## Removing poor fits completely via a threshold (good approach?)
#    wVec[wVec < 0.01] <- 0
    
    lenfl <- length(fctList)
    lenrl <- length(respLev)
    edEst <- matrix(NA, lenfl + 1, lenrl)
    edSe <- matrix(NA, lenfl + 1, lenrl)    
    
    ## Defining 'interval' argument for ED
    if (identical(interval, "kang"))
    {
        interval2 <- "delta"
    } else {
        interval2 <- "none"    
    }
    
    ## Calculating estimated ED values
    edMat <- ED(object, respLev, interval2, display = FALSE)
    edEst[1, ] <- as.vector((edMat)[, 1])
    edSe[1, ] <- as.vector((edMat)[, 2])
    
    if (identical(interval2, "delta"))
    {
        edCll <- matrix(NA, lenfl + 1, lenrl)
        edClu <- matrix(NA, lenfl + 1, lenrl)
        
        edCll[1, ] <- as.vector((edMat)[, 3])
        edClu[1, ] <- as.vector((edMat)[, 4])        
    }
    for (i in 1:lenfl)
    {
        edMati <- try(ED(update(object, fct = fctList[[i]]), respLev, interval2, display = FALSE), silent = TRUE)
        if (inherits(edMati, "try-error"))
        {
            edMati <- matrix(NA, length(respLev), 4)
        }
         
        edEst[i + 1, ] <- as.vector((edMati)[, 1])
        edSe[i + 1, ] <- as.vector((edMati)[, 2])
        if (identical(interval2, "delta"))
        {
            edCll[i + 1, ] <- as.vector((edMati)[, 3])
            edClu[i + 1, ] <- as.vector((edMati)[, 4])
        }
    }

    edVec <- apply(edEst * wVec, 2, sum, na.rm = na.rm)
    if (identical(interval, "none"))
    {
        retMat <- as.matrix(cbind(edVec))
        colnames(retMat) <- colnames(edMat)[1]
    }
    if (identical(interval, "buckland"))
    {
        seVec <- apply(sqrt(edSe^2 + (t(t(edEst) - apply(edEst, 2, mean, na.rm = na.rm)))^2) * wVec, 2, 
        sum, na.rm = na.rm)
### Thresholding                
#        iVec <- wVec < 0.01        
#        seVec <- apply(sqrt(edSe[iVec, ]^2 + (t(t(edEst[iVec, ]) - apply(edEst[iVec, ], 2, 
#        mean, na.rm = na.rm)))^2) * wVec[iVec], 2, sum, na.rm = na.rm)
        quantVal <- qnorm(1 - (1 - level)/2) * seVec
        retMat <- as.matrix(cbind(edVec, seVec, edVec - quantVal, edVec + quantVal))
        colnames(retMat) <- c(colnames(edMat)[c(1, 2)], "Lower", "Upper")
    }
    if (identical(interval, "kang"))
    {  
        retMat <- as.matrix(cbind(apply(edEst * wVec, 2, sum, na.rm = na.rm), 
        apply(edCll * wVec, 2, sum, na.rm = na.rm), 
        apply(edClu * wVec, 2, sum, na.rm = na.rm)))
        colnames(retMat) <- colnames(edMat)[c(1,3,4)]
    }   
    rownames(retMat) <- rownames(edMat)

    ## Constructing matrix of fit summaries 
    disMat <- as.matrix(cbind(edEst, wVec))
    colnames(disMat) <- c(paste("EC", rownames(edMat), sep = ""), "Weight")
    rownames(disMat) <- rownames(msMat)    
    if (display)
    {
        print(disMat)
        cat("\n") 
    }
    
#    resPrint(resMat, "Estimated effective doses", interval, "Model-averaging", display)
    if (extended)
    {
        return(list(estimates = retMat, fits = disMat))
    } else {
        retMat
    }
    
}