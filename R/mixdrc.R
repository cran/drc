"mixdrc" <- 
function(object, random, data, startVal)
{

#    if (missing(lambda)) {lambda <- 1}
    if (!is.null(object$"boxcox")) 
    {
        lambda <- object$"boxcox"[1]
        bcAdd <- object$"boxcox"[4]
    } else {
        lambda <- 1 
        bcAdd <- 1
    }

    ## Defining dose and response
    respVar <- ((eval(object$call[[2]][[2]], envir = data) + bcAdd)^lambda - 1)/lambda  # name of response variable
    doseVar <- eval(object$call[[2]][[3]], envir = data)
#    dataSet <- cbind(data, doseVar, respVar)
    dataSet <- data.frame(data, doseVar = doseVar, respVar = respVar)

    ## Getting parameter names
    parNames <- object$fct$names
    lenPN <- length(parNames)


    ## Defining dose-response function 
    if (inherits(object, "llogistic"))
    {
        if (lenPN == 3)
        {
            opfct <- function(DOSE,b,d,e)
            {
                ( (0 + (d-0)/(1+exp(b*(log(DOSE)-log(e)))) + bcAdd)^lambda - 1)/lambda
            }
            NLMEfor <- formula(respVar ~ opfct(doseVar, b, d, e))
        }
        if (lenPN == 4)
        {
            opfct <- function(DOSE,b,c,d,e)
            {
                ( (c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))) + bcAdd)^lambda - 1)/lambda
            }
            NLMEfor <- formula(respVar ~ opfct(doseVar, b, c, d, e))
        }
        if (lenPN < 5)
        {
            strVal <- "c + (d-c)/((1+exp(b*(log(DOSE)-log(e))))^f)"
            namesVec <- c("b", "c", "d", "e", "f")

            fctType <- deparse(object$call$fct[[1]])
            if (fctType %in% c("l3", "LL.3")) {fixedVec <- c(NA, 0, NA, NA, 1); indVec <- c(1, 3, 4)}
            if (fctType %in% c("l4", "LL.4")) {fixedVec <- c(NA, NA, NA, NA, 1); indVec <- c(1, 2, 3, 4)}

            fv <- eval(object$call$fct$fixed)
#            print(fv)
            if (!is.null(fv)) {fixedVec[indVec] <- fv}
#            print(fixedVec)

            retRC <- repChar(strVal, namesVec, fixedVec, keep=c("DOSE", "exp"))
            opfctStr <- retRC[[1]]
#            print(opfctStr)
            opfct <- eval(parse(text = opfctStr))
            
            NLMEfor <- eval(parse(text = retRC[[2]]))
#            print(NLMEfor)
        }
        
        
        if (lenPN == 5) {stop("Not implemented for l5")}       
    }
    assign("opfct", opfct, env = .GlobalEnv)  # assigning object to global environment    
#    print(opfct)
#   print(opfct(1, 1, 0.1, 1))



        ## Constructing list for fixed argument
        colEntry <- object$call$collapse

        fixedList <- list()
        if ( (!is.null(colEntry)) && (as.character(colEntry[[1]])=="list") )
        {
            for (i in 2:(1+lenPN))
            {
                fixedList[[i-1]] <- eval(parse(text=paste(parNames[i-1], deparse(object$call$collapse[[i]]), sep="")), envir = dataSet)
            }
        }
        if ( (!is.null(colEntry)) && (as.character(colEntry[[1]])=="data.frame"))
        {
            for (i in 2:(1+lenPN))
            {
                startStr <- as.character(object$call$collapse[[i]])
                if (length(grep("1", startStr)) == 0) 
                {
                    fili <- paste(parNames[i-1], paste("factor(", startStr, ")", sep=""), sep="~")    
                    fili <- paste(fili, "-1", sep="")
                } else {
                    fili <- paste(parNames[i-1], startStr, sep="~")    
                }
                
                fixedList[[i-1]] <- eval(parse(text=fili), envir = dataSet)
            }
        }
        if (is.null(colEntry)) 
        {
            for (i in 2:(1+lenPN))
            {
                if (!is.null(object$call$curve))
                {
                    facStr <- paste("factor(", deparse(object$call$curve), ")", sep="")            
                    fili <- paste(parNames[i-1], facStr, sep="~")
                    fili <- paste(fili, "-1", sep="")
                } else {
                    facStr <- "1"            
                    fili <- paste(parNames[i-1], facStr, sep="~")
                }
                
                fixedList[[i-1]] <- eval(parse(text=fili), envir = dataSet)
            }
        }

        
#        print(fixedList)


        ## Constructing list for random argument
        randomList <- eval(parse(text=random), envir = dataSet)
#        print(randomList)


        ## Searching for start values that yield convergence
        require(nlme, quietly = TRUE)
        
#        if (lenPN==3)
#        {
#            logist3 <- function(DOSE,b,d,e)
#            {
#                ( (0 + (d-0)/(1+exp(b*(log(DOSE)-log(e)))) + bcAdd)^lambda - 1)/lambda
#            }
#            assign("logist3", logist3, env=.GlobalEnv)  # assigning object to global environment
#
#            found <- FALSE
#            for (i in 1:30)
#            {
#
#                startVal <- c(i/10,coef(object)[-c(1)])
#
#                options(warn=-1)
#                modelNLME <- try(nlme(respVar~logist3(doseVar,b,d,e),
#                                      fixed = fixedList,
#                                      random = randomList,
#                                      start=startVal, na.action=na.omit, data=dataSet), silent=TRUE)
#
#                if (!inherits(modelNLME,"try-error")) {found <- TRUE; break}  # print(c(i)); stop("got it")}
#                options(warn=0)
#            }
#            rm(logist3, envir=.GlobalEnv)  # removing object from global environment            
#        }
#        if (lenPN==4)
#        {
#            logist4 <- function(DOSE,b,c,d,e)
#            {
#                ( (c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))) + bcAdd)^lambda - 1)/lambda
#            }
#            assign("logist4", logist4, env=.GlobalEnv)  # assigning object to global environment
#
#            found <- FALSE
#            for (i in 1:30)
#            {
#
#                startVal <- c(i/10,coef(object)[-c(1)])
#
#                options(warn=-1)
#                modelNLME <- try(nlme(respVar~logist4(doseVar,b,c,d,e),
#                                      fixed = fixedList,
#                                      random = randomList,
#                                      start=startVal, na.action=na.omit, data=dataSet), silent=TRUE)
#
#                if (!inherits(modelNLME,"try-error")) {found <- TRUE; break}  # print(c(i)); stop("got it")}
#                options(warn=0)
#            }
#            rm(logist4, envir=.GlobalEnv)  # removing object from global environment            
#        }
##        if (lenPN==5) {stop("Not implemented for l5")}       


    found <- FALSE
    
    if (missing(startVal))
    {
        coefNames <- names(coef(object))
        
#        print(NLMEfor)
#        print(fixedList)
#        print(randomList)
#        print((dataSet))
        
#        return(list(NLMEfor, fixedList, randomList, dataSet)) 
        
        for (i in 1:30)
        {
            startVal <- c(i/10, coef(object)[-c(1)])
            names(startVal) <- coefNames

            options(warn = -1)

#            print(dataSet)
#            print(fixedList)
#            print(randomList)
#            print(startVal)
            
#            modelNLME <- nlme(SLOPE ~ c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))),
#                             fixed = list(b~factor(HERBICIDE)-1, c~1, d~1, e~factor(HERBICIDE)-1),
#                             random = d~1|CURVE,
#                             start = as.vector(startVal), data = PestSci)
#            stop()
#            

#print(as.vector(startVal))

#            modelNLME <- nlme(SLOPE ~ c + (d-c)/(1+exp(b*(log(DOSE)-log(e)))),
#            fixed = list(b~factor(HERBICIDE)-1, c~1, d~1, e~factor(HERBICIDE)-1),
#            random = d~1|CURVE, start = as.vector(startVal), data = PestSci)
#
#stop()

            modelNLME <- try(nlme(NLMEfor,
                             fixed = fixedList,
                             random = randomList,
                             start = startVal, na.action = na.omit, data = dataSet), silent = TRUE)

            if (!inherits(modelNLME, "try-error")) 
            {
                found <- TRUE 
                break
            }
            options(warn = 0)
        }
    } else {
        print(fixedList)
        print(randomList)
    
        modelNLME <- try(nlme(NLMEfor,
                              fixed = fixedList,
                              random = randomList,
                              start = startVal, na.action = na.omit, data = dataSet)) ## , silent = TRUE)

        if (!inherits(modelNLME, "try-error")) {found <- TRUE}    
    }
    rm(opfct, envir = .GlobalEnv)  # removing object from global environment            
    
    if (!found) {stop("No convergence. The model may be too general.")}

#    }


    modelNLME$nlme <- modelNLME  # storing the original 'nlme' object
    modelNLME$data <- object$data
    modelNLME$df.residual <- df.residual(object)  # correct ???
    modelNLME$fct <- object$fct
    modelNLME$transformation <- object$transformation
    modelNLME$parmMat <- t(object$"pmFct"( as.vector((summary(modelNLME$nlme))$tTable[,1]) ))
    modelNLME$curve <- list( object$"pfFct"(t(modelNLME$parmMat)), object$curve[[2]])
#    modelNLME$parmMat <- t(modelNLME$parmMat)
    modelNLME$class <- "mixed logistic"
    modelNLME$parNames <- mdrcPNsplit(rownames(summary(modelNLME)$tTable), ".")
    modelNLME$base <- object

    class(modelNLME) <- c("mixdrc", "drc", class(modelNLME))
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
