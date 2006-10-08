"l5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), useDer = FALSE, w = FALSE)
{
    return(logistic(fixed = fixed, names = names, useDer = useDer, w = w))
}
