"anova.drclist" <- function(object, ...,  test = NULL)
{
    objects <- list(object, ...)
    if (length(objects) > 2) {stop("Only two models can be compared")}

    if (inherits(object, "bindrc"))  # the argument 'test="F"' is not used
    {
        obj1 <- objects[[1]]$loglik
        obj2 <- objects[[2]]$loglik
        rowNames <- c("Model 1", "Model 2")


        loglik <- c(obj1[3], obj2[3])
        dfModel <- c(obj1[4], obj2[4])
        testStat <- (2*abs(loglik[1] - loglik[2]))
        dfDiff <- c(NA, abs(dfModel[1] - dfModel[2]))

        pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
        testStat <- c(NA, testStat)

        headName <- "Analysis of deviance table\n"
        colNames <- c("ModelDf", "Loglik", "Df", "LR value", "p value")

    } else {

    ## Testing two models against each other
    obj1 <- objects[[1]]
    obj2 <- objects[[2]]
    rowNames <- c("1st model", "2nd model")
 
    sumObj1 <- summary(obj1)
    sumObj2 <- summary(obj2)

    if ( !(obj1$"type"==obj2$"type") ) {stop("The two models are based on different types on data")}
    if (obj1$"type" == "binomial" && (is.null(test)) ) {test <- "Chisq"}
    if (obj1$"type" == "continuous"  && (is.null(test)) ) {test <- "F"}

    if (!(test == "F"))  # chis-square based test
    {
        loglik <- c(logLik(obj1), logLik(obj2))  # c(sumObj1[[4]][1], sumObj2[[4]][1])
        dfModel <- c(attr(logLik(obj1), "df"),  attr(logLik(obj2), "df"))  # c(sumObj1[[4]][2], sumObj2[[4]][2])
        testStat <- (2*abs(loglik[1] - loglik[2]))
        dfDiff <- c(NA, abs(dfModel[1] - dfModel[2]))

        pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
        testStat <- c(NA, testStat)

        headName <- "ANOVA-like table\n"
        colNames <- c("ModelDf", "Loglik", "Df", "LR value", "p value")
        
    } else {  # F-test
        
        sumVec1 <- obj1[[4]]
        sumVec2 <- obj2[[4]]

        if (sumVec2[6]>sumVec1[6]) {sumTemp <- sumVec1; sumVec1 <- sumVec2; sumVec2 <- sumTemp; rowNames <- rowNames[c(2,1)]}

        loglik <- c(sumVec1[5],sumVec2[5])
        dfModel <- c(sumVec1[6], sumVec2[6])
#        loglikTemp <- c(sumVec1[5], sumVec2[5])
#        loglik <- c((loglikTemp[1]-loglikTemp[2])/(dfModel[1]-dfModel[2]), loglikTemp[2]/dfModel[2])
#        dfModel <- c(sumVec1[6], sumVec2[6])
#        dfDiff <- c((loglik[1]-loglik[2])/(dfModel[1]-dfModel[2]), loglik[2]/dfModel[2])
        dfDiff <- c(NA, dfModel[1]-dfModel[2])
#        testStat <- dfDiff[1]/dfDiff[2]
        testStat <- ((loglik[1]-loglik[2])/(dfModel[1]-dfModel[2]))/(loglik[2]/dfModel[2])

        pVal <- c(NA,1-pf(testStat, dfModel[1]-dfModel[2], dfModel[2]))
        testStat <- c(NA,testStat)

        headName <- "ANOVA table\n"
        colNames <- c("ModelDf", "RSS", "Df", "F value", "p value")
    }
#    dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)


    ## Specifying the models
    cat("\n")

    collapse1 <- obj1[[8]]$collapse
#    if (is.null(collapse1)) {collapse1 <- paste(deparse(obj1[[8]]$assayNo), "(for all parameters)")} else {collapse1 <- deparse(obj1[[8]]$collapse)}
    if (is.null(collapse1)) 
    {
        if (is.null(obj1[[8]]$curve))
        {
            collapse1 <- "1 (for all parameters)"
        } else {
            collapse1 <- paste(deparse(obj1[[8]]$curve), "(for all parameters)")
        }
    } else {
        collapse1 <- paste(deparse(obj1[[8]]$collapse), collapse="")
        collapse1 <- gsub("    ", "", collapse1, fixed=TRUE)  # removing extra spaces
    }

#    pos <- 1
#    if (is.data.frame(eval(collapse1))) 
    pos <- regexpr("data.frame(", collapse1, fixed=TRUE)
    if (pos>0)     
    {
#        collapse1 <- deparse(obj1[[8]]$collapse)
#        collapse1 <- substring(collapse1, pos+11, nchar(collapse1)-1)
        collapse1 <- substring(collapse1, 12, nchar(collapse1)-1)
    }
#    if (is.list(eval(collapse1)))    
    pos <- regexpr("list(", collapse1, fixed=TRUE)
    if (pos>0) 
    {
#        collapse1 <- deparse(obj1[[8]]$collapse)
#        collapse1 <- substring(collapse1, pos+5, nchar(collapse1)-1)
        collapse1 <- substring(collapse1, 6, nchar(collapse1)-1)
    }      

    collapse2 <- obj2[[8]]$collapse
#    if (is.null(collapse2)) {collapse2 <- paste(deparse(obj2[[8]]$assayNo), "(for all parameters)")} else {collapse2 <- deparse(obj2[[8]]$collapse)}
    if (is.null(collapse2)) 
    {
        if (is.null(obj2[[8]]$curve))
        {
            collapse2 <- "1 (for all parameters)"
        } else {
            collapse2 <- paste(deparse(obj2[[8]]$curve), "(for all parameters)")
        }
    } else {
        collapse2 <- paste(deparse(obj2[[8]]$collapse), collapse="")
        collapse2 <- gsub("    ", "", collapse2, fixed=TRUE)  # removing extra spaces
    }
    
#    if (is.data.frame(eval(collapse2)))     
    pos <- regexpr("data.frame(", collapse2, fixed=TRUE)    
    if (pos>0)
    {
#        collapse2 <- deparse(obj2[[8]]$collapse)
#        collapse2 <- substring(collapse2, pos+11, nchar(collapse2)-1)
        collapse2 <- substring(collapse2, 12, nchar(collapse2)-1)
    }

#    if (is.list(eval(collapse2))) 
    pos <- regexpr("list(", collapse2, fixed=TRUE)
    if (pos>0)
    {
#        collapse2 <- deparse(obj2[[8]]$collapse)
#        collapse2 <- substring(collapse2, pos+5, nchar(collapse2)-1)
        collapse2 <- substring(collapse2, 6, nchar(collapse2)-1)
    }               
    
    cat("1st model\n")
    cat(paste(" fct:      ", deparse(obj1[[8]]$fct), "\n", sep=""))
    cat(paste(" collapse: ", collapse1, "\n", sep=""))
    cat("2nd model\n")
    cat(paste(" fct:      ", deparse(obj2[[8]]$fct), "\n", sep=""))
    cat(paste(" collapse: ", collapse2, "\n\n", sep=""))
    }

    dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
    dimnames(dataFra) <- list(rowNames, colNames)
    structure(dataFra, heading = headName, class = c("anova", "data.frame"))
}
