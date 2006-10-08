"l3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), useDer = FALSE, w = FALSE)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed[1], 0, fixed[2:3], 1), names=c(names[1],"c",names[2:3],"f"), useDer=useDer, w=w ) )
}
