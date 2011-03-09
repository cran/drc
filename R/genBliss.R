"genBliss" <- function(
fixed = rep(NA, 8), names = c("b1", "b2", "c", "d", "e1", "e2", "f1", "f2"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 8
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
   
    ## Defining the nonlinear model function    
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        bVal <- parmMat[, 4]
        mVal <- parmMat[, 3]

#        LL5 <- function(doseVec, pMat)
#        {
#            LL.5()$fct(doseVec, pMat)
#        }
        mVal + (bVal - mVal) * (((mVal - LL.5()$fct(dose[ ,1], parmMat[, c(1, 3, 4, 5, 7), drop = FALSE])) / (mVal - bVal)) * ((mVal - LL.5()$fct(dose[ ,2], parmMat[, c(2, 3, 4, 6, 8), drop = FALSE])) / (mVal - bVal)))
    }    

    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initVal1 <- llogistic()$ssfct(dframe[, c(1, 3)])
            initVal2 <- llogistic()$ssfct(dframe[, c(2, 3)])
#            initVal <- c((llogistic()$ssfct(dframe))[c(1, 1:4, 4)], 1, 1) * c(-1, -1, rep(1, 6))
            initVal <- c(initVal1, initVal2)[c(1, 6, 2, 8, 4, 9, 5, 10)] 
    
            return(initVal[notFixed])
        }        
    }    
    
    ## Defining names
    names <- names[notFixed]

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL

    ## Returning list with model function information
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct = edfct, sifct = sifct,
    name = "genBliss",
    text = "Generalized Bliss", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "genBliss"
    invisible(returnList)
}
