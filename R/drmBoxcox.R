"drmBoxcox" <- function(boxcox, anovaFormula, dset)
{
    lambda <- 0

    isNumeric <- is.numeric(boxcox)
    if ( (isNumeric) || (is.logical(boxcox) && boxcox)  ) 
    {
        if (!isNumeric)
        {
            profLik <- boxcox(anovaFormula, lambda = seq(-2.6, 2.6, 1/10), plotit = FALSE, data = dset)  
            # boxcox in MASS
                
            maxIndex <- which.max(profLik$y)
            lambda <- (profLik$x)[maxIndex]
            boxcoxci <- drmBoxcoxCI(profLik)
        }
        if (isNumeric)
        {
            lambda <- boxcox
            boxcox <- TRUE
            boxcoxci <- c(NA, NA)                
        }
    } else {
        lambda <- NA
        boxcoxci <- c(NA, NA)
    }
    return(list(lambda, boxcoxci, boxcox))
}
