"g4" <-
function(names=c("b", "c", "d", "e"), useDer=FALSE)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==4)) {stop("Not correct 'names' argument")}
    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    return(gompertz(names=c(names), useDer=useDer, fixed=c(NA, NA, NA, NA)))
}
