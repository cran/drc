"l4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), useDer = FALSE)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( logistic(fixed = c(fixed, 1), names = c(names,"f"), useDer = useDer ) )
}
