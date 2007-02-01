"scaleX" <- 
function(dose, lowTol)
{
    if (max(dose) < lowTol)
    {
        scaleFact <- floor(log(mean(dose), 10))
         
        scaleFct <- function(parm, down = TRUE) 
                    {
                        if (down) {parm <- parm / (10^scaleFact)} else {parm <- parm * (10^scaleFact)}
                    }
    } else {
        scaleFct <- function(x, down = TRUE) {return(x)}  # the identity
    }            
    return(scaleFct)
}

"scaleX" <- function(scaleX)
{
    scaleFct <- function(parm, down = TRUE)
    {
        if (down) {parm/scaleX} else {parm*scaleX}
    }
    return(scaleFct)
}

"scaleY" <- function(scaleY)
{
    scaleFct <- function(parm, down = TRUE)
    {
        if (down) {parm/scaleY} else {parm*scaleY}
    }
    return(scaleFct)
}
