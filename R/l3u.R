"l3u" <-
function(fixed = c(NA, NA, NA), names = c("b", "c", "e"), useDer = FALSE, w = FALSE)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed[1:2], 1, fixed[3], 1), names=c(names[1:2],"d",names[3],"f"), useDer=useDer, w = w ) )
}
