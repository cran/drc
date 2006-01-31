"w4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), useDer = FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}

    return(weibull(fixed = fixed, names = names, useDer = useDer))
}
