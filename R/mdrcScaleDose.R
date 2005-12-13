"mdrcScaleDose" <- function(dose, lowTol = 1e-2)
{
    if (max(dose) < lowTol)
    {
        scaleFact <- floor(log(mean(dose), 10))
         
        scaleFct <- function(parm, down = TRUE) 
                    {
                        if (down) {parm <- parm / (10^scaleFact)} else {parm <- parm * (10^scaleFact)}
                    }
    } else {
        scaleFct <- function(x, down = TRUE) {return(x)}
    }            
    return(scaleFct)
}
