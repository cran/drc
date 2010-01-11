"EDhelper" <- function(parmVec, respl, reference, type, cond = TRUE)
{
    ## Converting absolute to relative
    if (type == "absolute") 
    {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
    } else {  
        p <- respl
    }
    ## Swapping p for an increasing curve
    if (cond)
    {
        if ((type == "relative") && (parmVec[1] < 0) && (reference == "control"))
        {
            p <- 100 - p
        }
    } else {
        if ((type == "relative") && (reference == "control"))
        {
            p <- 100 - p
        }
    }
    p
}