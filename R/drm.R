"drm" <-
function(formula, curve, pmodels, weights, data = NULL, subset, fct, 
adjust = c("none", "bc1", "bc2", "vp"), bc = NULL, bcAdd = 0, 
type = c("continuous", "binomial", "Poisson", "survival"),
start, na.action = na.fail, hetvar = NULL, robust = "mean", logDose = NULL, 
fctList = NULL, control = drmc(), lowerl = NULL, upperl = NULL)
{
    ## Matching 'adjust' argument
    adjust <- match.arg(adjust)
    type <- match.arg(type)

    
    ## Loading MASS
    require(MASS, quietly = TRUE)  # used for boxcox and ginv


    ## Setting na.action option
    options(na.action = deparse(substitute(na.action)))


    ## Setting control parameters
    useD <- control$"useD"
    constrained <- control$"constr"
    maxDose <- control$"maxDose"
    maxIt <- control$"maxIt"
    optMethod <- control$"method"
    relTol <- control$"relTol"
    warnVal <- control$"warnVal"
    zeroTol <- control$"zeroTol"
    bcConstant <- bcAdd  # the bcAdd argument overrules the setting via the control argument
    rmNA <- control$"rmNA"
    errorMessage <- control$"errorm"
    noMessage <- control$"noMessage"
    
    
    ## Setting warnings policy
    options(warn = warnVal)


    ## Setting adjustment
    if (adjust == "none") {boxcox <- FALSE; varPower <- FALSE}
    if (adjust == "bc1") {boxcox <- TRUE; varPower <- FALSE}
    if (adjust == "vp") {boxcox <- FALSE; varPower <- TRUE}
    if ( (!is.null(bc)) && (is.numeric(bc))) {boxcox <- bc}
        
    if (!(robust == "mean"))
    {
        boxcox <- FALSE
        varPower <- FALSE
    }


    ## Handling 'start' argument
    if (missing(start)) {selfStart <- TRUE} else {selfStart <- FALSE}


    ## Handling 'fct' argument
    if ( (!is.list(fct)) && (!is.function(fct)) ) {stop("No function or list given in argument 'fct'")}
    if (is.function(fct)) 
    {
        fct <- fct2list(fct, 2)
    }
    
    ## Converting a user specified list
    if (is.null(names(fct))) {fct$"fct" <- fct[[1]]; fct$"ssfct" <- fct[[2]]; fct$"names" <- fct[[3]]}
    
    if (!is.function(fct$"fct")) 
    {
        stop("First entry in list to 'fct' NOT a function")
    } else {
        drcFct <- fct$"fct"
    }
    
    if (is.null(fct$"ssfct")) {noSSfct <- TRUE} else {noSSfct <- FALSE}
    if ((!is.function(fct$"ssfct")) && selfStart)
    {
        stop("Neither self starter function nor starting values provided")
    } else {
        ssfct <- fct$"ssfct"
    }
    
    if (is.null(fct$"names") || (!is.character(fct$"names"))) 
    {
        stop("Parameter names (as vector a strings) are NOT supplied")
    } else {
        parNames <- fct$"names" 
        numNames <- length(parNames)
    }
    
    ## Coercing two arguments in 'ssfct' into one argument
    lenASS <- length(formals(ssfct))
    if (lenASS > 1)
    {
        ssTemp <- ssfct
        ssfct <- function(dataset) {ssTemp(dataset[, head(1:lenASS, -1)], dataset[, lenASS])}
    }

    ## Checking whether or not first derivates are supplied    
    if ( (useD) && (is.function(fct$"deriv1")) )
    {
        dfct1 <- fct$"deriv1"  # deriv1  # [[4]]
#        drcDer2 <- fct$deriv2  # [[5]]
    } else {
        dfct1 <- NULL
    }

    ## Checking whether or not second derivates are supplied    
    if ( (useD) && (is.function(fct$"deriv2")) )
    {
        dfct2 <- fct$"deriv2"
    } else {
        dfct2 <- NULL
    }
    
    fct$"anovaYes"$"bin" <- NULL
    fct$"anovaYes"$"cont" <- TRUE


    ## Handling the 'formula', 'curve' and 'data' arguments
    anName <- deparse(substitute(curve))  # storing name for later use
    if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)

    mf <- match.call(expand.dots = FALSE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curve", "data", "subset", "na.action", "weights", "hetvar"), nmf, 0) 

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())  #, globalenv())
    mt <- attr(mf, "terms")
        
    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
#    print(dose)
    resp <- model.response(mf, "numeric")
    origResp <- resp  # in case of transformation of the response    
    lenData <- length(resp)
    xDim <- ncol(as.matrix(dose))
    dimData <- xDim + 1  # dimension of dose plus 1 dimensional response
    
    varNames <- names(mf)
    varNames <- varNames[c(2:dimData,1)]


    ## Retrieving weights
    weights <- model.weights(mf)
    if (is.null(weights))
    {
        weights <- rep(1, lenData)
    }


    ## Extracting variable for heterogeneous variances
    vvar <- model.extract(mf, "hetvar")   
    
    
    ## Finding indices for missing values
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}


    ## Handling situation where no 'curve' or 'pmodels' argument is given
    assayNo <- model.extract(mf, "curve")   
    pmodelsList <- list()

    if (is.null(assayNo)) {assayNo <- rep(1,lenData)}
    uniqueNames <- unique(assayNo)
    colOrder <- order(uniqueNames)
    uniqueNames <- as.character(uniqueNames)
    
    if (missing(pmodels)) 
    {
        pmodels <- as.data.frame(matrix(assayNo, lenData, numNames))
        
        if (length(unique(assayNo)) == 1) 
        {
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- matrix(1, lenData, 1)
        }} else {
            modelMat <- model.matrix(~factor(assayNo)-1)  # no intercept term
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- modelMat
                pmodelsList[[i]] <- pmodelsList[[i]][, colOrder]
        }}   
    } else {

    ## Handling a list or data.frame argument of 'pmodels'
    if (is.null(data)) 
    {
        pmodels <- eval(substitute(pmodels), envir = .GlobalEnv)
    } else {
        pmodels <- eval(substitute(pmodels), envir = data, enclos = parent.frame())
    }
   
    if (is.data.frame(pmodels))
    {
        lenCol <- ncol(pmodels)
        pmodelsMat <- matrix(0, lenData, lenCol)    
    
        for (i in 1:lenCol) 
        {
            if (length(unique(pmodels[,i])) == 1) 
            {
                pmodelsList[[i]] <- matrix(1, lenData, 1)
                pmodelsMat[,i] <- rep(1, lenData)    
            }
            else {
                mf <- eval(model.frame(~factor(pmodels[,i]) - 1), parent.frame())  # converting to factors
                mt <- attr(mf, "terms")    
    
                mf2 <- model.matrix(mt, mf)
                ncmf2 <- ncol(mf2)

                mf3 <- removeMI(mf2)
                pmodelsList[[i]] <- mf3
                pmodelsMat[, i] <- mf3%*%c(1:ncmf2)
            }
        }
           
    } else {

        if (is.list(pmodels))
        {   
            lenCol <- length(pmodels)
            pmodelsMat <- matrix(0, length(resp), lenCol)
    
            for (i in 1:lenCol) 
            {
                if (paste(as.character(pmodels[[i]]), collapse = "") == "~1") 
                {
                    pmodelsList[[i]] <- matrix(1, lenData, 1)
                    pmodelsMat[,i] <- rep(1, lenData)
                } else
                {
                    mf <- eval(model.frame(pmodels[[i]], data=data), parent.frame())   
                    mt <- attr(mf, "terms")    
                        
                    mf2 <- model.matrix(mt, mf)
                    ncmf2 <- ncol(mf2)

                    mf3 <- removeMI(mf2)                    
                    pmodelsList[[i]] <- mf3  
                    
                    pmodelsMat[,i] <- mf3%*%c(1:ncmf2)                    
                }
            }
        }
    }     
    pmodelsOld <- pmodels
    pmodels <- as.data.frame(pmodelsMat)  # pmodelsMat not used any more
    }
    

    ## Re-setting na.action
    options(na.action = "na.omit")  # the default


    ## Transforming dose value if they are provided as log dose
    if ( !is.null(logDose) && is.numeric(logDose) ) 
    {
       origDose <- dose
       dose <- logDose^dose
    }


    ## Re-enumerating the levels in 'assayNo' and 'pmodels'
    assayNoOld <- assayNo


    ## Defining helper function     
    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    }

    assayNo <- colConvert(assayNoOld)  # curves enumerated from 1 to numAss
    assayNames <- as.character(unique(assayNoOld))
    numAss <- length(assayNames)
    for (i in 1:numNames) {pmodels[, i] <- colConvert(pmodels[, i])}


    ## Handling one-dimensional x     
    if (xDim == 1)
    {
        ## Detecting control measurements 
        uniqueDose <- lapply(tapply(dose, assayNoOld, unique), length)
        udNames <- names(uniqueDose[uniqueDose == 1])
        if (length(udNames) > 0) 
        {
            cm <- udNames
            if (!noMessage) {cat(paste("Control measurements detected for level: ", udNames, "\n", sep = ""))}
        } else {
            cm <- NULL
        }


        ## Defining ANOVA model
        bcc <- rep(bcConstant, lenData)    
        if (numAss > 1) 
        {
            anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor*factor(assayNo)
            alternative <- 2
        } else {
            anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor
            alternative <- 1
        }
         

        ## Checking whether there is enough df to perform Box-Cox transformation
        if ( boxcox && ( (lenData - numAss*length(unique(dose))) < lenData/10) )
        {
            if (boxcox) {warning("Box-Cox transformation based on clustering of dose values", call. = FALSE)}
            doseFactor <- factor(cutree(hclust(dist(dose), method = "average"), lenData/3))  
            # constructing groups containing roughly 3 observations 


            ## Re-defining ANOVA model
            if (numAss > 1) 
            {
                anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor*factor(assayNo)
                dset <- data.frame(doseFactor, resp, assayNo, bcc)
                alternative <- 2
            } else {
                anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor
                dset <- data.frame(doseFactor, resp, bcc)
                alternative <- 1
            }
        } else {
            doseFactor <- factor(dose)
        }
        dset <- data.frame(dose, doseFactor, resp, assayNo, bcc)

         
        ## Fitting ANOVA model
        if (type == "continuous")
        {
            testList <- drmLOFls()
            if (varPower) {testList <- drmLOFvp()}
            if (!is.null(vvar)) {testList <- drmLOFhv()}
        }
        if (type == "binomial")
        {
            testList <- drmLOFbinomial()
        }
        if (type == "Poisson")
        {
            testList <- drmLOFPoisson()
        }

        
#        if (varPower) {testList <- mdrcVp(anovaYes = TRUE)} else {testList <- drmEMls(anovaYes = TRUE)}
#        if (!is.null(vvar)) {testList <- mdrcHetVar(anovaYes = TRUE)}

        gofTest <- testList$"gofTest"            
        lofTest <- testList$"anovaTest"
        if (!is.null(lofTest))
        {
            anovaModel0 <- lofTest(anovaFormula, dset)
        } else {
            anovaModel0 <- NULL
            alternative <- 0
        }


        ## Applying the Box-Cox transformation (lambda is defined here!)
        bcResult <- drmBoxcox(boxcox, anovaFormula, dset)  
        lambda <- bcResult[[1]]
        boxcoxci <- bcResult[[2]] 
        boxcox <- bcResult[[3]]      

#        lambda <- 0
#        isNumeric <- is.numeric(boxcox)
#        if ( (isNumeric) || (is.logical(boxcox) && boxcox)  ) 
#        {
#            if (!isNumeric)
#            {
#                profLik <- boxcox(anovaFormula, lambda = seq(-2.6, 2.6, 1/10), plotit = FALSE, data = dset)  
#                # boxcox in MASS
#                
#                maxIndex <- which.max(profLik$y)
#                lambda <- (profLik$x)[maxIndex]
#                boxcoxci <- drmBoxcoxCI(profLik)
#            }
#            if (isNumeric)
#            {
#                lambda <- boxcox
#                boxcoxci <- c(NA, NA)                
#            }
#        } else {
#            lambda <- NA
#            boxcoxci <- c(NA, NA)
#        }
 

        ## Using self starter 
        if (!noSSfct)
        {
            ## Calculating initial estimates for the parameters using the self starter
            startMat <- matrix(0, numAss, numNames)
            doseresp <- data.frame(dose, origResp)
    
            for (i in 1:numAss)
            {
                indexT1 <- (assayNo == i)
                if (any(indexT1)) 
                {
                    isfi <- is.finite(dose)  # removing infinite dose values
                    logVec <- indexT1 & isfi
                    startMat[i, ] <- ssfct(doseresp[logVec, ])  # ssfct(dose[logVec], origResp[logVec] )
                } else {
                    startMat[i, ] <- rep(NA, numNames)
                }

                ## Identifying a dose response curve only consisting of control measurements
                if (sum(!is.na(startMat[i, ])) == 1) {upperPos <- (1:numNames)[!is.na(startMat[i, ])]}
            }
#            colMat <- matrix(0, numNames, numAss)
#            maxParm <- rep(0, numNames)  # storing the max number of parameters
        }

        
    ## Handling multi-dimensional x   
    } else {  
        alternative <- NULL
        anovaModel0 <- NULL
#        anovaModel <- NULL
        gofTest <- NULL
        
        if (is.numeric(boxcox))
        {
            lambda <- boxcox                  
            boxcoxci <- c(NA, NA)
        } else {
            lambda <- NA                  
            boxcoxci <- NULL
        }
        

        ## Using self starter
        if (!noSSfct)
        {
            ## Calculating initial estimates for the parameters using the self starter
            startMat <- matrix(0, numAss, numNames)
            doseresp <- data.frame(dose, origResp)
    
            for (i in 1:numAss)
            {
                indexT1 <- (assayNo == i)
                if (any(indexT1)) 
                {
                    startMat[i, ] <- ssfct(doseresp[indexT1, ])  # ssfct(dose[indexT1], origResp[indexT1])
                } else {
                    startMat[i, ] <- rep(NA, numNames)
                }
                
                ## Identifying a dose response curve only consisting of control measurements
                if (sum(!is.na(startMat[i,]))==1) {upperPos <- (1:numNames)[!is.na(startMat[i,])]}
            }
#            colMat <- matrix(0, numNames, numAss)
#            maxParm <- rep(0, numNames)  # storing the max number of parameters
        }
    }


    ## Finding parameters for the control measurements which will not be estimated
    pmodelsList2 <- list()
    for (i in 1:numNames)
    { 
        colNames <- colnames(pmodelsList[[i]])

        if ( (!is.null(cm)) && (!is.null(colNames)) ) 
        {
            accm <- as.character(cm)
            pos <- grep(accm, colNames)
            if (length(pos) == 0) 
            {
                candCol <- pmodelsList[[i]][, 1]
                if ( !(length(assayNoOld[candCol==1])==0) && (all(assayNoOld[candCol==1] == accm)) )
                {
                    pos <- 1  # the control measurements correspond to the "Intercept" term
                }
            }  
        } else {pos <- numeric(0)}


        ## Defining 'pmodelsList2' from 'pmodelsList'
        if ((length(pos) > 0) && !(upperPos == i) )
        {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]][, -pos])  # column is removed
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])  # column is kept
        } 
    }  
   
    
    
    ## Constructing vectors 'ncclVec' and 'parmPos' used below
    ncclVec <- rep(0, numNames)
    for (i in 1:numNames)
    {
        ncclVec[i] <- ncol(pmodelsList2[[i]])  # ncol(as.matrix(pmodelsList2[[i]]))
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])


    ## Constructing vector of initial parameter values
    startVecList <- list()
    if(!noSSfct)
    {
        nrsm <- nrow(startMat)
        for (i in 1:numNames)
        {
            sv <- rep(0, max(nrsm, ncclVec[i]))
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]
            
            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)
            
            startVecList[[i]] <- sv[indVec]
        }
        startVec <- unlist(startVecList)        
    } else {
        startVec <- start  # no checking if no self starter function is provided!!!
    }


    ## Checking the number of start values provided
    if (!selfStart && !noSSfct) 
    {
        lenReq <- length(startVec)
        if (length(start) == lenReq) 
        {
            startVec <- start
        } else {
            stop(paste("Wrong number of initial parameter values. ", lenReq, " values should be supplied", sep=""))
        }
    }


    ## Converting parameters
    if (selfStart)
    {
        startVec <- drmConvertParm(startVec, startMat, assayNo, pmodelsList2) 
    }


    ## Constructing parameter names
    pnList <- drmParNames(numNames, parNames, pmodelsList2)
    parmVec <- pnList[[1]]
    parmVecA <- pnList[[2]]
    parmVecB <- pnList[[3]]
    

    ## Defining function which converts parameter vector to parameter matrix            
    parmMatrix <- matrix(0, lenData, numNames)
    parm2mat <- function(parm)
    {
#        parmMatrix <- matrix(0, lenData, numNames)
#print(parm)
        for (i in 1:numNames)
        {
#           print(as.matrix(pmodelsList2[[i]]))
           parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmPos[i] + 1:ncclVec[i]]
        }
        return(parmMatrix)
    }
        

    ## Defining non-linear function
    if (!is.null(fctList))
    {
        ivList <- list()
        ivList2 <- list()        
        matList <- list()
        svList <- list()
        for (i in 1:numAss)
        {
            indexT1 <- (assayNo == i)
            isfi <- is.finite(dose)  # removing infinite dose values

            ivList[[i]] <- indexT1
#            svList[[i]] <- fctList[[i]]$"ssfct"( doseresp[(indexT1 & isfi), ] )
            logVec <- indexT1 & isfi
            svList[[i]] <- fctList[[i]]$"ssfct"(doseresp[logVec, ])  # dose[logVec], origResp[logVec])
            matList[[i]] <- c( sum(indexT1), length(svList[[i]]) )
            
            ivList2[[i]] <- match(fctList[[i]]$names, fct$names)             
        }

    
        posVec <- rep(0, numAss)
        for (i in 1:numAss)
        {
            posVec[i] <- matList[[i]][2]
        }
        posVec <- cumsum(posVec)
        posVec <- c(0, posVec)
#        print(posVec)
    
        drcFct1 <- function(dose, parm)
        {
            retVec <- rep(0, lenData)
            for (i in 1:numAss)
            {
                iVec <- ivList[[i]]
#                print(matList[[i]])
                pMat <- matrix(parm[(posVec[i]+1):posVec[i+1]], matList[[i]][1], matList[[i]][2], byrow = TRUE) 
#                print(pMat)
#                print(length(dose[iVec]))
#                print(dim(pMat))
                retVec[iVec] <- fctList[[i]]$"fct"( dose[iVec], pMat )
            }
            return(retVec)
        }
        
        startVec <- as.vector(unlist(svList))
    } else {

        drcFct1 <- function(dose, parm)
        {
            drcFct(dose, parm2mat(parm))
        }
    }

    if (is.null(cm)) 
    {
        multCurves <- function(dose, parm)
        {          
           drcFct1(dose, parm)  # fctList           
        }
    } else {
        iv <- assayNoOld == cm
        niv <- !iv
        fctEval <- rep(0, lenData)

        multCurves<-function(dose, parm) 
        {
            parmVal <- parm2mat(parm)           
            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , drop = FALSE])

            fctEval
        }
    }


    ## Defining objective function
    robustFct <- drmRobust(robust, match.call(), lenData, length(startVec))
    

    ## Box-Cox transformation is applied
    if (boxcox)
    {
#        varPower <- FALSE  # not both boxcox and varPower at the same time

        ## Defining Box-Cox transformation function
        bcfct <- function(x, lambda, bctol, add = bcConstant)
        {
            if (abs(lambda) > bctol)
            {
                return(((x+add)^lambda-1)/lambda)
            } else {
                return(log(x+add))    
            }
        }
        
        ## Setting the tolerance for Box-Cox transformation being the logarithm transformation 
        ##  (same as in boxcox.default in MASS package)
        bcTol <- 0.02 
        
        resp <- bcfct(resp, lambda, bcTol)

        multCurves2 <- function(dose, parm)
        {
            bcfct(multCurves(dose, parm), lambda, bcTol)
        } 
    } else {multCurves2 <- multCurves}
    
    
    ## Defining first derivative (if available)
    if (!is.null(dfct1))
    {
        dmatfct <- function(dose, parm)
        {
            dfct1(dose, parm2mat(parm))
        }
    } else {
        dmatfct <-NULL
    }
    
   
    ## Defining estimation method
    scaleXConstant <- 1
    scaleYConstant <- 1
    
    if (type == "continuous")
    {
        ## Ordinary least squares estimation
        estMethod <- drmEMls(dose, resp, multCurves2, startVec, robustFct, weights, rmNA, dmf = dmatfct, 
        scaleX = scaleXConstant, scaleY = scaleYConstant)
        
        if (varPower)
        {        
            estMethod <- drmEMvp(dose, resp, multCurves2)  # mdrcVp(dose, resp, multCurves2)
            lenStartVec <- length(startVec)
            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
            parmVec <- c(parmVec, "Sigma", "Power")

        }
                  
        if (!is.null(vvar))
        {
            estMethod <- mdrcHetVar(dose, resp, multCurves2, vvar)
            lenStartVec <- length(startVec)
            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
            parmVec <- c(parmVec, as.character(unique(vvar)))        
        }
    }
    if (type == "binomial")
    {
        estMethod <- drmEMbinomial(dose, resp, multCurves2, startVec, robustFct, weights, rmNA)    
    
    }
    if (type == "Poisson")
    {
        estMethod <- drmEMPoisson(dose, resp, multCurves2, startVec)
    }       
    opfct <- estMethod$opfct            


    ## Re-fitting the ANOVA model to incorporate Box-Cox transformation (if necessary)
    if (!is.na(lambda))
    {
        dset <- data.frame(dose, doseFactor, resp, assayNo, bcc)  # dataset with new resp values        
        anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)            
#        anovaModel <- anovaModel0$"anovaFit"
    }


    ## Defining lower and upper limits of parameters
    if (constrained)
    {
        if (!is.null(lowerl)) 
        {
            if (!is.numeric(lowerl) || !((length(lowerl) == sum(ncclVec)) || (length(lowerl) == numNames)))
            {
                stop("Not correct 'lowerl' argument")
            } else {
                if (length(lowerl) == numNames) 
                {
                    lowerLimits <- rep(lowerl, ncclVec)
                } else {
                    lowerLimits <- lowerl
                }
            } 
        }

        if (!is.null(upperl)) 
        {
            if (!is.numeric(upperl) || !((length(upperl) == sum(ncclVec)) || (length(upperl) == numNames)))
            {
                stop("Not correct 'upperl' argument")
            } else {
                if (length(upperl) == numNames) 
                {
                    upperLimits <- rep(upperl, ncclVec)
                } else {
                    upperLimits <- upperl
                }
            } 
        }
        
        if (all(!is.finite(lowerLimits)) && all(!is.finite(upperLimits))) 
        {
            stop("No constraints are imposed via 'lowerl' and 'upperl' arguments")
        }
    }


    ## Optimising
    
    ## Setting derivatives
    opdfctTemp <- estMethod$"opdfct1"
    appFct <- function(x, y){tapply(x, y, sum)}   
    
    if (!is.null(opdfctTemp))
    {
        opdfct1 <- function(parm)
        {
#            print(as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo)))
            as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo))
        }
    } else {
        opdfct1 <- NULL
    }  


    ## Manipulating before optimisation
#    print(startVec)  # 1
        
    ## Scaling x values
    sxInd <- fct$"sxInd"
    sxYN <- !is.null(sxInd) && ((max(dose)<1e-2) || (min(dose)>1e2) || (diff(range(dose))>1e2) )
    if ( sxYN && (is.null(fctList)) )
    {
#        if (!is.null(fctList))
#        {
#            parmIndX <- rep(0, numAss)
#            for (i in 1:numAss)
#            {
#                parmIndX[i] <- fctList[[i]]$"sxInd"
#            }
#            parmIndX <- cumsum(parmIndX)
#        } else {
            parmIndX <- parmPos[sxInd] + 1:ncclVec[sxInd]        
#        }
    
        scaleXConstant <- median(dose)
        sxFct <- scaleX(scaleXConstant)  # , scaleX(dose, maxDose)        
#        dose <- sxFct(dose)
        startVec[parmIndX] <- sxFct(startVec[parmIndX])
    }
#    print(startVec)  # 2

    ## Scaling y values
    syInd <- fct$"syInd"
    lensy <- length(syInd)
    parmIndY <- list()
    syYN <- !is.null(syInd) && ((max(resp)<1e-2) || (min(resp)>1e2) || (diff(range(resp))>1e2)) 
    if ( syYN && (is.null(fctList)) )
    {
#        if (!is.null(fctList))
#        {
#            parmIndY <- rep(0, numAss)
#            for (i in 1:numAss)
#            {
#                parmIndY[[i]] <- fctList[[i]]$"syInd"
#            }
#            parmIndY <- cumsum(as.vector(unlist(parmIndY)))
#        } else {
            for (i in 1:lensy)
            {
                parmIndY[[i]] <- parmPos[syInd[i]] + c(1:ncclVec[syInd[i]])
            }
            tempPIY <- as.vector(unlist(parmIndY))
            parmIndY <- tempPIY
#        }
        scaleYConstant <- median(resp) 
        syFct <- scaleY(scaleYConstant)
        startVec[parmIndY] <- syFct(startVec[parmIndY])
    }
#    print(startVec)  # 3


    ## Testing nonlinear function
#    print(startVec)
#    print(multCurves2(dose, startVec))
#    print(opfct(startVec))
#    print(dose)
#    print(resp)
    

    ## Optimising
    nlsFit <- drmOpt(opfct, opdfct1, startVec, optMethod, constrained, warnVal, 
    upperLimits, lowerLimits, errorMessage, maxIt, relTol, parmVec = parmVec) 
    
    nlsFit$ovalue <- nlsFit$value    


    ## Manipulating after optimisation
    
    ## Adjusting for scaling of y values
    if ( syYN && (is.null(fctList)) )
    {
        nlsFit$value <- syFct(syFct(nlsFit$value, down = FALSE), down = FALSE)
        startVec[parmIndY] <- syFct(startVec[parmIndY], down = FALSE)
        nlsFit$par[parmIndY] <- syFct(nlsFit$par[parmIndY], down = FALSE)

        scaleFct1 <- function(hessian) 
                     {
                         newHessian <- hessian
                         newHessian[, parmIndY] <- syFct(newHessian[, parmIndY], down = FALSE)
                         newHessian[parmIndY, ] <- syFct(newHessian[parmIndY, ], down = FALSE)
                         return(newHessian)
                     }                
    } else {
        scaleFct1 <- function(x) {x}    
    }

    
    ## Adjusting for scaling of x values
    if ( sxYN && (is.null(fctList)) )  # (!is.null(sxInd))
    {
#        dose <- sxFct(dose, down = FALSE)
        startVec[parmIndX] <- sxFct(startVec[parmIndX], down = FALSE)
        nlsFit$par[parmIndX] <- sxFct(nlsFit$par[parmIndX], down = FALSE)

        scaleFct2 <- function(hessian) 
                     {
                         newHessian <- scaleFct1(hessian)
                         newHessian[, parmIndX] <- sxFct(newHessian[, parmIndX], down = FALSE)
                         newHessian[parmIndX, ] <- sxFct(newHessian[parmIndX, ], down = FALSE)
                         return(newHessian)
                     }                
    } else {
        scaleFct2 <- function(hessian) 
        {
            scaleFct1(hessian)
        }
    }


    ## Handling variance parameters
    varParm <- NULL
        
    if (varPower)
    {
        varParm <- list(type = "varPower", index = 1:lenStartVec)        
    }
    if (!is.null(vvar))
    {
        varParm <- list(type = "hetvar", index = 1:lenStartVec)
    }


    # Testing against the ANOVA (F-test)
    nlsSS <- nlsFit$value
    nlsDF <- lenData - length(startVec)


    ## Constructing a plot function
        
    ## Picking parameter estimates for each curve. Does only work for factors not changing within a curve!
    if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames==cm)]} else {iVec <- 1:numAss}
    
    pickCurve <- rep(0, length(iVec))
    for (i in iVec)
    {
       pickCurve[i] <- (1:lenData)[assayNo==i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)

    fixedParm <- (estMethod$"parmfct")(nlsFit)
#    print(nlsFit$par)
#    print(fixedParm)
    parmMat[iVec,] <- parm2mat(fixedParm)[pickCurve, ]

    if(!is.null(fctList))
    {
         parmMat <- matrix(NA, numAss, numNames)
         for (i in 1:numAss)
         {
             parmMat[i, ivList2[[i]]] <- fixedParm[(posVec[i]+1):posVec[i+1]]
         }
    }    
        
    if (!is.null(cm))
    {
#        conPos <- upperPos
#        print(conPos)
        parmMat[-iVec, upperPos] <- parm2mat(fixedParm)[assayNoOld==cm, , drop=FALSE][1, upperPos]  
        # 1: simply picking the first row
    }
    rownames(parmMat) <- assayNames


    pmFct <- function(fixedParm)
    {
        if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames == cm)]} else {iVec <- 1:numAss}
    
        if (!is.null(cm))
        {
#            conPos <- conList$"pos"
            parmMat[-iVec, upperPos] <- parm2mat(fixedParm)[assayNoOld==cm, , drop=FALSE][1, upperPos]  
            # 1: simply picking the first row
        }
        rownames(parmMat) <- assayNames
    
        return(parmMat)
    }
    parmMat <- pmFct( (estMethod$"parmfct")(nlsFit) )


    ## Constructing design matrix allowing calculations for each curve
    colPos <- 1
    rowPos <- 1
#    Xmat <- matrix(0, numAss*numNames, length(nlsFit$par))
#    Xmat <- matrix(0, numAss*numNames, length(fixedParm))


#    if (!is.null(fctList)) {omitList <- list()}
#    for (i in 1:numNames)
#    {
#        indVec <- iVec
#        lenIV <- length(indVec)
#
#        nccl <- ncol(pmodelsList2[[i]])  # min(maxParm[i], ncol(pmodelsList2[[i]]))        
#
#        XmatPart <- matrix(0, lenIV, nccl)
#        k <- 1
#        if (!is.null(fctList)) {omitVec <- rep(TRUE, lenIV)}
#        for (j in indVec)
#        {
#            if (!is.null(fctList))
#            {
#                parPresent <- !is.na(match(i, ivList2[[j]]))
#                omitVec[k] <- parPresent
#            }
#            
#            XmatPart[k, ] <- (pmodelsList2[[i]])[(1:lenData)[assayNo == j][1], 1:nccl]
#            k <- k + 1
#        }
#        if (!is.null(fctList))
#        {
#            XmatPart <- XmatPart[omitVec, , drop = FALSE]
#            nccl <- nccl - sum(!omitVec)
#            omitList[[i]] <- omitVec
#        }
#
#        Xmat[rowPos:(rowPos+lenIV-1), colPos:(colPos+nccl-1)] <- XmatPart
#        colPos <- colPos + nccl
#        rowPos <- rowPos + lenIV
#    }
#    Xmat <- Xmat[1:(rowPos-1), 1:(colPos-1)]


    ## Defining the plot function    
    pfFct <- function(parmMat)
    {
        plotFct <- function(dose)
        {
            if (xDim == 1) {lenPts <- length(dose)} else {lenPts <- nrow(dose)}

            curvePts <- matrix(NA, lenPts, numAss)
            for (i in 1:numAss)
            {
                if (!is.null(fctList)) 
                {
                    drcFct <- fctList[[i]]$"fct"
                    numNames <- matList[[i]][2]    
                }
            
                if (i%in%iVec)
                {
#                    parmChosen <- parmMat[i, ]
                    parmChosen <- parmMat[i, complete.cases(parmMat[i, ])]  # removing NAs 
#                    print(parmChosen)                   
                    
                    parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow = TRUE)
#                    print(parmMat2)
                    curvePts[, i] <- drcFct(dose, parmMat2)
                } else { curvePts[, i] <- rep(NA, lenPts)}
            }
            return(curvePts)
        }
    
        return(plotFct)    
    }
#    print(parmMat)
    plotFct <- pfFct(parmMat)
#    plotFct(0:10)


    ## Computation of fitted values and residuals
    predVec <- multCurves2(dose, fixedParm)    
    resVec <- resp - predVec
    resVec[is.nan(predVec)] <- 0    


    diagMat <- matrix(c(predVec, resVec), lenData, 2)
    colnames(diagMat) <- c("Predicted values", "Residuals")


    ## Adjusting for robust estimation: MAD based on residuals, centered at 0, is used as scale estimate    
    if (robust%in%c("median", "trimmed", "tukey", "winsor"))
    {
        nlsFit$value <- (mad(resVec, 0)^2)*nlsDF 
    }
#    if (robust=="winsor")
#    {
#        K <- 1 + length(startVec)*var(psi.huber(resVec/s, deriv=1))
#    }
    if (robust%in%c("lms", "lts"))  # p. 202 i Rousseeuw and Leroy: Robust Regression and Outlier Detection
    {  
        scaleEst <- 1.4826*(1+5/(lenData-length(nlsFit$par)))*sqrt(median(resVec^2))                                 
        w <- (resVec/scaleEst < 2.5)
        nlsFit$value <- sum(w*resVec^2)/(sum(w)-length(nlsFit$par))    
    }

    
    ## Adding meaningful names for robust methods
    robust <- switch(robust, median="median", trimmed="metric trimming", tukey="Tukey's biweight", 
                             winsor="metric Winsorizing", lms="least median of squares",
                             lts="least trimmed squares")


    ## Collecting summary output
    sumVec <- c(lambda, NA, NA, NA, nlsSS, nlsDF, lenData, alternative)
    sumList <- list(lenData = lenData, alternative = alternative, df.residual = lenData - length(startVec))


    ## The function call
    callDetail <- match.call()
    if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}


    ## The data set
    if (!is.null(logDose)) 
    {
        dose <- origDose
    }
    dataSet <- data.frame(dose, origResp, assayNo, assayNoOld, weights)
    names(dataSet) <- c(varNames, anName, anName, "weights")


    ## Box-Cox information
    bcVec <- c(lambda, boxcoxci)
    if (all(is.na(bcVec))) {bcVec <- NULL}
    if (!is.null(bcVec)) {bcVec <- c(bcVec, bcAdd)}


    ## Evaluating goodness-of-fit test
    if (!is.null(gofTest)) {gofTest <- gofTest(resp, weights, predVec, sumList$"df.residual")}


    ## Adjusting in case 'fctList' is specified
    if (!is.null(fctList))
    {
        omitAllVec <- as.vector(unlist(omitList))

        parmVec <- parmVec[omitAllVec] 
        parmVecA <- parmVecA[omitAllVec] 
        parmVecB <- parmVecB[omitAllVec] 
        
        orderVec <- match(as.vector(parmMat), nlsFit$par)
        orderVec <- orderVec[complete.cases(orderVec)]       
        
        nlsFit$par <- nlsFit$par[orderVec]
        nlsFit$hessian <- nlsFit$hessian[orderVec, orderVec]
    }


    ## Constructing an index matrix for use in ED and SI
    hfct1 <- function(x)  # helper function
    {
        uniVec <- unique(x[!is.na(x)])
        rv <- rep(NA, length(x))
        for (i in 1:length(uniVec))
        {
            rv[abs(x-uniVec[i]) < 1e-12] <- i
        }
        rv
    }
    hfct2 <- function(x)
    {
        length(unique(x))
    }
#    parmMat <- t(parmMat)
    mat1 <- t(apply(t(parmMat), 1, hfct1))  # , 1:ncol(parmMat)))
    cnccl <- head(cumsum(ncclVec), -1)
#    mat2 <- mat1
    if (nrow(mat1) == 1) {mat1 <- t(mat1)}  # in case of only one curve
    mat1[-1, ] <- mat1[-1, ] + cnccl

    

    ## Returning the fit
    returnList <- list(varParm, nlsFit, list(plotFct, logDose), sumVec, startVec, list(parmVec, parmVecA, parmVecB), 
    diagMat, callDetail, dataSet, t(parmMat), fct, robust, bcVec, estMethod, lenData-length(startVec), 
    anovaModel0, gofTest, sumList, scaleFct2, pmFct, pfFct, type, mat1, logDose)
    names(returnList) <- c("varParm", "fit", "curve", "summary", "start", "parNames", "predres", "call", "data", 
    "parmMat", "fct", "robust", "boxcox", "estMethod", "df.residual", "anova", "gofTest", 
    "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat", "logDose")
    class(returnList) <- c("drc", class(fct))

    return(returnList)
}
