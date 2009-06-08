"EDhelper" <- function(parmVec, respl, reference, type)
{
    ## Converting absolute to relative
    if (type == "absolute") 
    {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
    } else {  
        p <- respl
    }
    ## Swapping p for increasing curve
    if ( (type == "relative") && (parmVec[1] < 0) && (reference == "control") )
    {
        p <- 100 - p
    }
    p
}