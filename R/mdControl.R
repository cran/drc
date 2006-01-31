"mdControl" <- function(bcAdd = 0, constr = FALSE, errorm = TRUE, maxIt = 500, method = "BFGS", noMessage = FALSE,
                        relTol = 1e-7, rmNA = FALSE, warnVal = -1, zeroTol = 0)
{
    return(list(bcAdd = bcAdd,
                constr = constr,
                errorm = errorm, 
                maxIt = maxIt, 
                method = method,
                noMessage = noMessage,
                relTol = relTol,
                rmNA = rmNA, 
                warnVal = warnVal, 
                zeroTol = zeroTol))
}
