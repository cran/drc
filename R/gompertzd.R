"gompertzd" <- function(
lowerc = c(-Inf, -Inf), upperc = c(Inf, Inf), fixed = c(NA, NA), names = c("a", "b"))
{   
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        innerT1 <- parmMat[, 2]*dose
        innerT2 <- parmMat[, 1]/parmMat[, 2]*(exp(innerT1) - 1)
        parmMat[, 1]*exp(innerT1 - innerT2)
    }


    ## Defining the self starter function
    ssfct <- function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]    
    
        aVal <- max(y)
        bVal <- 1
        
        return(c(aVal, bVal))
    }
   
    ## Defining names
    names <- names[notFixed]
        
    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        help1 <- fct(dose, parm)
        help2 <- exp(parmMat[, 2]*dose)
        help3 <- help2 - 1
        help4 <- parmMat[, 1]/parmMat[, 2]
        
        deriva <- help1*(1/parmMat[, 1] - help3/parmMat[, 2])
        derivb <- help1*(dose + help4*help3/parmMat[, 2] + help4*help2*dose)
        
        cbind(deriva, derivb)[, notFixed]       
    }
        
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        fct(x, parm)*(parmMat[, 2] - parmMat[, 1]*exp(parmMat[, 2]*x))
    }

    ## Setting the limits
    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

    ## Defining the ED function
    edfct <- NULL
if (FALSE)
{    
    function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
        if (type == "absolute") 
        {
            p <- 100*((parmVec[3] - respl)/(parmVec[3] - parmVec[2]))
        } else {  
            p <- respl
        }
        if ( (parmVec[1] < 0) && (reference == "control") )
        {
            p <- 100 - p
        }
    
        tempVal <- log((100-p)/100)
        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])

        EDder <- 
        EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
        0, 0, 1/parmVec[4], 
        exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))

        return(list(EDp, EDder[notFixed]))
    }
}

    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, lowerc = lowerLimits, upperc = upperLimits)
    
    class(returnList) <- "gompertzd"
    invisible(returnList)
}

if (FALSE)
{
"LL.2" <-
function(fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, 1, fixed[2], 1), names = c(names[1], "c", "d", names[2], "f"), ...) )
}

l2 <- LL.2

"LL.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f"), ...) )
}

l3 <- LL.3

"LL.3u" <-
function(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1:2], 1, fixed[3], 1), names = c(names[1:2], "d", names[3], "f"), ...) )
}

l3u <- LL.3u

"LL.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)  #, reps = FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed, 1), names = c(names, "f"), ...) )  # , reps = reps) )
}

l4 <- LL.4

"LL.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return(llogistic(fixed = fixed, names = names, ...))
}

l5 <- LL.5
}
