"drmc" <- 
function(constr = FALSE, errorm = TRUE, maxDose = 1e-1, maxIt = 500, method = "BFGS", 
noMessage = FALSE, relTol = 1e-7, rmNA = FALSE, useD = FALSE, trace = FALSE, warnVal = -1, zeroTol = 0)
{
    return(list(
                constr = constr,
                errorm = errorm,
                maxDose = maxDose,
                maxIt = maxIt, 
                method = method,
                noMessage = noMessage,
                relTol = relTol,
                rmNA = rmNA, 
                useD = useD,
                trace = trace,
                warnVal = warnVal, 
                zeroTol = zeroTol))
}
