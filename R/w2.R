"w2" <-
function(fixed = c(NA, NA), names = c("b", "e"), useDer = FALSE)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull(fixed = c(fixed[1], 0, 1, fixed[2]), names = c(names[1], "c", "d", names[2]), useDer = useDer))
}
