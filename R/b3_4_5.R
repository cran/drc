"b3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( boltzmann(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f")) )
}


"b4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( boltzmann(fixed = c(fixed, 1), names = c(names, "f")) )
}


"b5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))
{
    return( boltzmann(fixed = fixed, names = names) )
}
