"mixdrc" <- 
function(object, random, lambda=1, data)
{

#    if (missing(lambda)) {lambda <- 1}


    ## Defining dose and response
    respVar <- (eval(object$call[[2]][[2]], envir=data)^lambda - 1)/lambda  # name of response variable
    doseVar <- eval(object$call[[2]][[3]], envir=data)
    dataSet <- cbind(data, doseVar, respVar)


    ## Getting parameter names
    parNames <- object$fct$names
    lenPN <- length(parNames)


    ## Defining dose-response function 
    if (inherits(object, "logistic"))
    {

#        if (lenPN==3)
#        {
#            logist3 <- function(DOSE,b,d,e)
#            {
#                ( (0 + (d-0)/(1+exp(b*(log(DOSE)-log(e)))) )^lambda - 1)/lambda
#            }
#            assign("logist3", logist3, env=.GlobalEnv)
#        }
#        if (lenPN==4)
#        {
#            logist4 <- function(DOSE,b,c,d,e)
#            {
#                ( (c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))) )^lambda - 1)/lambda
#            }
#            assign("logist4", logist4, env=.GlobalEnv)
#        }
#        if (lenPN==5) {stop("Not implemented for l5")}       


        ## Constructing list for fixed argument
        colEntry <- object$call$collapse

        fixedList <- list()
        if (as.character(colEntry[[1]])=="list")
        {
            for (i in 2:(1+lenPN))
            {
                fixedList[[i-1]] <- eval(parse(text=paste(parNames[i-1], deparse(object$call$collapse[[i]]), sep="")))
            }
        }
        if (as.character(colEntry[[1]])=="data.frame")
        {
            for (i in 2:(1+lenPN))
            {
                fili <- paste(parNames[i-1], as.character(object$call$collapse[[i]]), sep="~")
                if (length(grep("1", fili))==0) {fili <- paste(fili, "-1", sep="")}
                
                fixedList[[i-1]] <- eval(parse(text=fili))
            }
        }
        if (is.null(colEntry)) 
        {
            for (i in 2:(1+lenPN))
            {
                fili <- paste(parNames[i-1], as.character(object$call$curve), sep="~")
                fili <- paste(fili, "-1", sep="")
                
                fixedList[[i-1]] <- eval(parse(text=fili))
            }
        }

        
#        print(fixedList)


        ## Constructing list for random argument
        randomList <- eval(parse(text=random))
#        print(randomList)


        ## Searching for start values that yield convergence
        require(nlme, quietly = TRUE)
        
        if (lenPN==3)
        {
            logist3 <- function(DOSE,b,d,e)
            {
                ( (0 + (d-0)/(1+exp(b*(log(DOSE)-log(e)))) )^lambda - 1)/lambda
            }
            assign("logist3", logist3, env=.GlobalEnv)

            found <- FALSE
            for (i in 1:30)
            {

                startVal <- c(i/10,coef(object)[-c(1)])

                options(warn=-1)
                modelNLME <- try(nlme(respVar~logist3(doseVar,b,d,e),
                                      fixed = fixedList,
                                      random = randomList,
                                      start=startVal, na.action=na.omit, data=dataSet), silent=TRUE)

                if (!inherits(modelNLME,"try-error")) {found <- TRUE; break}  # print(c(i)); stop("got it")}
                options(warn=0)
            }
            rm(logist3, envir=.GlobalEnv)  # removing object from global environment            
        }
        if (lenPN==4)
        {
            logist4 <- function(DOSE,b,c,d,e)
            {
                ( (c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))) )^lambda - 1)/lambda
            }
            assign("logist4", logist4, env=.GlobalEnv)

            found <- FALSE
            for (i in 1:30)
            {

                startVal <- c(i/10,coef(object)[-c(1)])

                options(warn=-1)
                modelNLME <- try(nlme(respVar~logist4(doseVar,b,c,d,e),
                                      fixed = fixedList,
                                      random = randomList,
                                      start=startVal, na.action=na.omit, data=dataSet), silent=TRUE)

                if (!inherits(modelNLME,"try-error")) {found <- TRUE; break}  # print(c(i)); stop("got it")}
                options(warn=0)
            }
            rm(logist4, envir=.GlobalEnv)  # removing object from global environment            
        }
        if (lenPN==5) {stop("Not implemented for l5")}       
        
        if (!found) {stop("No convergence. The model may be too general.")}

    }



    modelNLME$class <- "logistic model with random effects"
    modelNLME$parNames <- mdrcPNsplit(rownames(summary(modelNLME)$tTable), ".")

    class(modelNLME) <- c("mixdrc", class(modelNLME))
    return(modelNLME)
}


#sumExtend <- function(object)
#{
#    cat("\n")
#    cat(paste("A '", object$class, "' was fit.\n", sep = ""))
#    cat("\n")
#    cat("Parameter estimates:\n\n")
#
#    resultMat <- as.matrix(summary(object)$tTable[,c(1,2,4,5)]) 
#    printCoefmat(resultMat)
#
#    varComp <- matrix(as.numeric(VarCorr(object)[,1]))
#    colnames(varComp) <- "Variance"
#    rownames(varComp) <- rownames(VarCorr(object))
#
#    cat("\nEstimated variance components:\n\n")
#    printCoefmat(varComp)
#
#
#    sumObj <- summary(object)
#
#    ll <- logLik(object)
#    loglik <- ll[1] 
#    degfre <- sumObj$dims$N - attr(ll, "df")
#
#    estimates <-  as.vector(sumObj$coefficients$fixed)
#    parNames <- rownames(resultMat)
#    varMat <- sumObj$varFix
#
#
#    ## Defining return list
#    retList <- list(varComp, varMat, resultMat, c(loglik, degfre), parNames)
#
#    names(retList) <- c("varComp", "varMat", "estimates", "loglik", "parNames") 
#    class(retList) <- c("summary.drc")
#    return(retList)
#}
