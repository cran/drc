"iceLoewe2.1" <- function (
fixed = rep(NA, 8), names = c("b1", "b2", "c1", "c2", "d", "e1", "e2", "f"), ssfct = NULL)
{
    numParm <- 8
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
        for (k in 1:25)
        {
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
        applyImplicitFct <- function(parmVec)
        {
            if ((!is.finite(parmVec[6])) && (!is.finite(parmVec[7])))
            {
                retVal <- parmVec[5]
            } else {
                implicitFct <- function(oVal)
                {
                    (1/oVal - 1)^(1/parmVec[1]) / parmVec[6] +
                    (1/oVal - 1)^(1/parmVec[2]) / parmVec[7] +
                    (1/oVal - 1)^(mean(1/parmVec[1:2])) * abs(parmVec[8]) / (parmVec[6] * parmVec[7]) - 1
                }
                bisection <- try(bisec(implicitFct, 0, 1), silent = TRUE)
                if (inherits(bisection, "try-error"))
                {
                  occuprate <- NA
                }
                else {
                  occuprate <- bisection$root
                }
                retVal <- (parmVec[3] + occuprate * (parmVec[5] - parmVec[3])) * (1 / occuprate - 1)^(1/parmVec[1]) / parmVec[6] +
                          (parmVec[4] + occuprate * (parmVec[5] - parmVec[4])) * (1 / occuprate - 1)^(1/parmVec[2]) / parmVec[7] +
                          (parmVec[5] - sign(parmVec[8]) * mean(c(abs(parmVec[5] - parmVec[3]), abs(parmVec[5] - parmVec[3]))) +
                              occuprate * parmVec[8] * mean(c(abs(parmVec[5] - parmVec[3]), abs(parmVec[5] - parmVec[4]))) *
                              (1 / occuprate - 1)^(mean(1/parmVec[1:2]))) / (parmVec[6] * parmVec[7])
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
            initval <- c((llogistic()$ssfct(dframe))[c(1, 1:2, 2:4, 4)], 0.5) * c(-1, -1, rep(1, 6))
            return(initval[notFixed])
        }
    }

    names <- names[notFixed]

    deriv1 <- NULL

    deriv2 <- NULL

    edfct <- NULL

    sifct <- NULL

    returnList <- list(fct = fct, ssfct = ssfct, names = names,
        deriv1 = deriv1, deriv2 = deriv2, edfct = edfct, sifct = sifct,
        name = "iceLoewe2.1",
        text = "iceLoewe with diff. maximal responses, single parameter",
        noParm = sum(is.na(fixed)))

    class(returnList) <- "iceLoewe2"

    invisible(returnList)
}

