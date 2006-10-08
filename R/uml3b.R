"uml3b" <-
function(names = c("b", "d", "e", "f"), useDer = FALSE)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}
    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    return(umlogistic(names = c(names[1], "c", names[2:4]), useDer = useDer, fixed = c(NA, 0, NA, NA, NA), alpha = 0.5))
}
