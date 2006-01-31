"contAnovaTest" <- function()
{
    anovaTest <- function(formula, ds)
    {
        anovaFit <- lm(formula, data = ds)
        if (df.residual(anovaFit) > 0)
        {
            return(list(test = "F", anovaFit = anovaFit))
        } else {
            return(NULL)
        }
    }
    return(anovaTest)
}
