"g3" <-
function(names=c("b", "d", "e"), useDer=FALSE)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==3)) {stop("Not correct 'names' argument")}
    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    return(gompertz(names=c(names[1], "c", names[2:3]), useDer=useDer, fixed=c(NA, 0, NA, NA)))
}
