"ml4c" <-
function(names=c("b", "c", "d", "e", "f"), useDer=FALSE)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==5)) {stop("Not correct 'names' argument")}
    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    return(mlogistic(names=names, useDer=useDer, fixed=c(NA, NA, NA, NA, NA), alpha=0.25))
}
