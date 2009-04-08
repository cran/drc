"drmLOFls" <- function()
{
    ## Defining lack-of-fit/goodness-of-fit tests
    anovaTest <- contAnovaTest()
    gofTest <- NULL

    return(list(anovaTest = anovaTest, gofTest = gofTest))
}