"iceLoewe.1" <- function (
fixed = rep(NA, 7), names = c("b1", "b2", "c", "d", "e1", "e2", "f"), ssfct = NULL)
{
    numParm <- 7
    if (!is.character(names) | !(length(names) == numParm))
    {
        stop("Not correct 'names' argument")
    }
    if (!(length(fixed) == numParm))
    {
        stop("Not correct 'fixed' argument")
    }
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    bisec <- function(fu, fuLow, fuHigh)
    {
        for (k in 1:25) {
            fuMiddle <- 0.5 * (fuLow + fuHigh)
            if (fu(fuMiddle) > 0) {
                fuHigh <- fuMiddle
            } else {
                fuLow <- fuMiddle
            }
        }
        list(root = fuMiddle)
    }
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[, 5] <- parmMat[, 5] / dose[, 1]
        parmMat[, 6] <- parmMat[, 6] / dose[, 2]

        applyImplicitFct <- function(parmVec)
        {
            if ((!is.finite(parmVec[5])) && (!is.finite(parmVec[6])))
            {
                retVal <- parmVec[4]
            } else {
                implicitFct <- function(oVal)
                {
                    (1/oVal - 1)^(1/parmVec[1]) / parmVec[5] +
                    (1/oVal - 1)^(1/parmVec[2]) / parmVec[6] +
                    ((1/oVal - 1)^((1/parmVec[1]+1/parmVec[2])/2) / (parmVec[5] * parmVec[6])) * (parmVec[7]) - 1
                }
                bisection <- try(bisec(implicitFct, 0, 1), silent = TRUE)
                if (inherits(bisection, "try-error"))
                {
                  occuprate <- NA
                }
                else {
                  occuprate <- bisection$root
                }
                retVal <- (parmVec[3] + occuprate * (parmVec[4] - parmVec[3])) *
                          ((1 / occuprate - 1)^(1/parmVec[1]) / parmVec[5] + (1 / occuprate - 1)^(1/parmVec[2]) / parmVec[6] +
#                          parmVec[7] * ((1 / occuprate - 1)^((1/parmVec[1]+1/parmVec[2])/2)) / (parmVec[5] * parmVec[6]))
                          parmVec[7] * ((1 / occuprate - 1)^((1/parmVec[1]+1/parmVec[2])/2)) / (parmVec[5] * parmVec[6]))
	          }
            return(retVal)
        }
        apply(parmMat, 1, applyImplicitFct)
    }

    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            startLL.d1 <- as.vector(coef(drm(dframe[, c(3,1)], fct = LL.4())))
            startLL.d2 <- as.vector(coef(drm(dframe[, c(3,2)], fct = LL.4())))

            if (startLL.d1[1] < 0)  # condition in terms of standard "drc" parameter "b"
            {
                initVal <- c(startLL.d1, startLL.d2, 0)[c(1, 5, 3, 2, 4, 8, 9)]
            } else {
                initVal <- c(startLL.d1, startLL.d2, 0)[c(1, 5, 2, 3, 4, 8, 9)] * c(-1, -1, rep(1, 5))
            }

            return(initVal[notFixed])
        }
    }
    
    names <- names[notFixed]
    
    deriv1 <- NULL
    
    deriv2 <- NULL
    
    edfct <- NULL
    
    sifct <- NULL

    returnList <- list(fct = fct, ssfct = ssfct, names = names,
        deriv1 = deriv1, deriv2 = deriv2, edfct = edfct, sifct = sifct,
        name = "iceLoewe.1new",
        text = "iceLoewe with common maximal response, two parameters",
        noParm = sum(is.na(fixed)))
        
    class(returnList) <- "iceLoewe.1" #?
    
    invisible(returnList)
}
