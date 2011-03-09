"genBliss2" <- function(
fixed = rep(NA, 9), names = c("b1", "b2", "c1", "c2", "d", "e1", "e1", "f1", "f2"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 9
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

        bVal <- parmMat[, 5]
        mVal <- pmax(parmMat[, 3], parmMat[, 4])

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
            initVal <- c(initVal1, initVal2)[c(1, 6, 2, 7, 8, 4, 9, 5, 10)] 
    
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
    name = "genBliss2",
    text = "Generalized Bliss2", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "genBliss2"
    invisible(returnList)
}
