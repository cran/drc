"multdrc" <-
function(formula, curve, collapse, weights, data = NULL, boxcox = FALSE, bcAdd = 0, varPower = FALSE, startVal, fct = l4(), 
         na.action = na.fail, hetvar = NULL, robust = "mean", type = "continuous", cm = NULL, logDose = NULL, 
         fctList = NULL, control = mdControl())
{
    require(MASS, quietly = TRUE)  # used for boxcox and ginv

    ## Setting na.action option
    options(na.action=deparse(substitute(na.action)))


    ## Setting control parameters
    constrained <- control$"constr"
    maxDose <- control$"maxDose"
    maxIt <- control$"maxIt"
    optMethod <- control$"method"
    relTol <- control$"relTol"
    warnVal <- control$"warnVal"
    zeroTol <- control$"zeroTol"
    bcConstant <- control$"bcAdd"
    bcConstant <- bcAdd  # the bcAdd argument overrules the setting via the control argument
    rmNA <- control$"rmNA"
    errorMessage <- control$"errorm"
    noMessage <- control$"noMessage"
    
    
    ## Setting warnings policy
    options(warn = warnVal)
    

#    design <- FALSE
    upperPos <- 0
#    zeroTol <- -0.1  # for use with separate control group; any number below 0 will work
#    zeroTol <- 5e-3

    ## Tolerance for Box-Cox transformation being the logarithm transformation (same as in boxcox.default in MASS package)
    bcTol <- 0.02  


    ## Setting distribution
    if (!is.null(fct$"dist")) {type <- fct$"dist"}
    if (type == "binomial") {boxcox <- FALSE; varPower <- FALSE}
    if (type == "continuous")
    {
        if (boxcox) {varPower <- FALSE}
        if (varPower) {boxcox <- FALSE}
        
        if (!(robust == "mean"))
        {
            boxcox <- FALSE
            varPower <- FALSE
        }
    }


    ## Handling 'startVal' argument
    if (missing(startVal)) {selfStart <- TRUE} else {selfStart <- FALSE}


    ## Handling 'fct' argument
    if (!is.list(fct)) {stop("No list given as argument to 'fct'")}
     
#    if (!is.function(fct[[1]])) {stop("First entry in list to 'fct' is NOT a function")} else {drcFct <- fct[[1]]}
    
    ## Converting a user specified list
    if (is.null(names(fct))) {fct$"fct" <- fct[[1]]; fct$"ssfct" <- fct[[2]]; fct$"names" <- fct[[3]]}
    
    if (!is.function(fct$"fct")) 
    {
        stop("First entry in list to 'fct' is NOT a function")
    } else {
        drcFct <- fct$"fct"
    }
    
#    if (is.null(fct[[2]])) {noSSfct <- TRUE} else {noSSfct <- FALSE}
    if (is.null(fct$"ssfct")) {noSSfct <- TRUE} else {noSSfct <- FALSE}
    if ((!is.function(fct$"ssfct")) && selfStart)
    {
        stop("No self starter function nor starting values are provided")
    } else {
        ssfct <- fct$"ssfct"
    }
    
#    if (!is.function(fct[[2]]) & selfStart) {stop("Second entry in list to 'fct' is NOT a function")} else {ssfct <- fct[[2]]}
#    if (!is.function(fct[[2]]) && selfStart && !noSSfct) {stop("Second entry in list to 'fct' is NOT a function")} else {ssfct <- fct[[2]]}   
#    if (!is.character(fct[[3]])) {stop("Third entry in list to 'fct' is NOT parameter names (strings)")} else {
    
    if (is.null(fct$"names") || (!is.character(fct$"names"))) 
    {
        stop("Parameter names (as vector a strings) are NOT supplied")
    } else {
        parNames <- fct$"names" 
        numNames <- length(parNames)
    }
    
    if (is.function(fct$deriv1) & is.function(fct$deriv2))
    {
        derFlag <- TRUE
        drcDer <- fct$deriv1  # [[4]]
        drcDer2 <- fct$deriv2  # [[5]]
    } else {
        derFlag <- FALSE
    }
    
    fct$"anovaYes"$"bin" <- NULL
    fct$"anovaYes"$"cont" <- TRUE


    ## Handling the 'formula', 'curve' and 'data' arguments
    anName <- deparse(substitute(curve))  # storing name for later use
    if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)

    mf <- match.call(expand.dots = FALSE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curve", "data", "na.action", "weights", "hetvar"), nmf, 0) 

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())
    mt <- attr(mf, "terms")
        
    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
#    print(dose)
    resp <- model.response(mf, "numeric")
    origResp <- resp  # in case of transformation of the response    
    lenData <- length(resp)
    doseDim <- ncol(as.matrix(dose))
    dimData <- doseDim + 1  # dimension of dose plus 1 dimensional response
    
    varNames <- names(mf)
    varNames <- varNames[c(2:dimData,1)]


    ## Retrieving weights
    weights <- model.weights(mf)
    if ( (is.null(weights)) & (type == "continuous") ) 
    {
        weights <- rep(1, lenData)
    }
    if ( (is.null(weights)) & (type == "binomial") ) 
    {
        stop("Argument 'weights' needs to be specified")
    }

    ## Extracting variable for heterogeneous variances
    vvar <- model.extract(mf, "hetvar")   
    
    
    ## Finding indices for missing values
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}


    ## Handling situation where no 'curve' or 'collapse' argument is given
    assayNo <- model.extract(mf, "curve")   
    collapseList <- list()

    if (is.null(assayNo)) {assayNo <- rep(1,lenData)}
    uniqueNames <- unique(assayNo)
    colOrder <- order(uniqueNames)
    uniqueNames <- as.character(uniqueNames)
    
    if (missing(collapse)) 
    {
        collapse <- as.data.frame(matrix(assayNo,lenData,numNames))
        
        if (length(unique(assayNo))==1) 
        {
            for (i in 1:numNames) 
            {
                collapseList[[i]] <- matrix(1, lenData, 1)
        }} else {
            modelMat <- model.matrix(~factor(assayNo)-1)  # no intercept term
            for (i in 1:numNames) 
            {
                collapseList[[i]] <- modelMat
                collapseList[[i]] <- collapseList[[i]][, colOrder]
        }}   
    } else {

    ## Handling a list or data.frame argument of 'collapse'
    if (is.null(data)) 
    {
        collapse <- eval(substitute(collapse), envir=.GlobalEnv)
    } else {
        collapse <- eval(substitute(collapse), envir=data, enclos=parent.frame())
    }
   
    if (is.data.frame(collapse))
    {
        lenCol <- ncol(collapse)
        collapseMat <- matrix(0, lenData, lenCol)    
    
        for (i in 1:lenCol) 
        {
            if (length(unique(collapse[,i])) == 1) 
            {
                collapseList[[i]] <- matrix(1, lenData, 1)
                collapseMat[,i] <- rep(1, lenData)    
            }
            else {
                mf <- eval(model.frame(~factor(collapse[,i])-1), parent.frame())  # converting to factors
#                mf <- eval(model.frame(~factor(collapse[,i])), parent.frame())  # converting to factors
                mt <- attr(mf, "terms")    
    
                mf2 <- model.matrix(mt, mf)
                ncmf2 <- ncol(mf2)

                mf3 <- removeMI(mf2)
                collapseList[[i]] <- mf3
#                collapseList[[i]] <- collapseList[[i]][, colOrder]                
                collapseMat[,i] <- mf3%*%c(1:ncmf2)
                
#                                    
#                collapseList[[i]] <- collapseList[[i]][, colOrder]
#                collapseList[[i]] <- collapseList[[i]]    
#                collapseMat[,i] <- model.matrix(mt, mf)%*%c(1:ncmf2)    
            }
        }   
    } else {

        if (is.list(collapse))
        {
#            design <- TRUE
    
            lenCol <- length(collapse)
            collapseMat <- matrix(0, length(resp), lenCol)
    
            for (i in 1:lenCol) 
            {
                if (paste(as.character(collapse[[i]]),collapse="")=="~1") 
                {
                    collapseList[[i]] <- matrix(1, lenData, 1)
                    collapseMat[,i] <- rep(1, lenData)
                } else
                {
                    mf <- eval(model.frame(collapse[[i]], data=data), parent.frame())   
                    mt <- attr(mf, "terms")    
                        
                    mf2 <- model.matrix(mt, mf)
#                    print(mf2$assign)
#                    assignCol <- as.logical(attr(mf2, "assign"))
#                    mf2 <- mf2[, assignCol]
                    ncmf2 <- ncol(mf2)

                    mf3 <- removeMI(mf2)
#                    print(assignCol)
#                    mf3 <- mf3[, assignCol]
                    
                    collapseList[[i]] <- mf3  
#                    print(dim(mf3))                  
                    
#                    attrTL <- attr(mt, "term.labels")
#                    print(attrTL)
#                    if ( (length(attrTL)==1) && (regexpr(":",attrTL)>0) && !(regexpr("-1",attrTL)>0) ) {collapseList[[i]] <- mf3[, -c(1), drop=FALSE]}
#                    print(dim(collapseList[[i]]))

#                    print(dim(collapseList[[i]]))
                    
                    collapseMat[,i] <- mf3%*%c(1:ncmf2)                    
#
#                    if (is.null(missingIndices))
#                    {
#                        collapseList[[i]] <- mf2
#                        collapseMat[,i] <- mf2%*%c(1:ncmf2)                    
#                    } else {
#                        collapseList[[i]] <- mf2[-missingIndices, ]
#                        collapseMat[,i] <- mf2[-missingIndices, ]%*%c(1:ncmf2)
#                    }
                }
            }
        }
    }     
    collapseOld <- collapse
    collapse <- as.data.frame(collapseMat)  # collapseMat not used any more
    }
    

    ## Re-setting na.action
    options(na.action = "na.omit")  # the default


    ## Transforming dose value if they are provided as log dose
    if (!is.null(logDose)) 
    {
       origDose <- dose
#       dose <- exp(dose)
       dose <- logDose^dose
    }


    ## Re-enumerating the levels in assayNo and collapse
    assayNoOld <- assayNo
#    numAss <- length(unique(assayNoOld))
    
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
    for (i in 1:numNames) {collapse[, i] <- colConvert(collapse[, i])}


    ## Handling one-dimensional dose data     
    if (doseDim == 1)
    {
        ## Detecting control measurements 
        uniqueDose <- lapply(tapply(dose, assayNoOld, unique), length)
        udNames <- names(uniqueDose[uniqueDose == 1])
        if ( (is.null(cm)) && (length(udNames) > 0) ) 
        {
            cm <- udNames
            if (!noMessage) {cat(paste("Control measurements detected for level: ",udNames, "\n", sep=""))}
        }


        ## Fitting the ANOVA model
        if (type == "binomial")
        {
            if (numAss > 1) 
            {
                anovaFormula <- cbind(resp*weights, weights - resp*weights) ~ factor(dose)*factor(assayNo)
                alternative <- 2
            } else {
                anovaFormula <- cbind(resp*weights, weights - resp*weights) ~ factor(dose)
                alternative <- 1
            }
            dset <- data.frame(dose, resp, weights, assayNo)

#            if (!is.null(fct$"anovaTest"))
            testList <- mdrcBinomial(anovaYes = TRUE)
            gofTest <- testList$"gofTest"                
                        
            if (!is.null(fct$"anovaYes"$"bin")) 
            {
#                anovaModel0 <- (fct$"anovaTest")(anovaFormula, dset)
                anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)
            } else {
                anovaModel0 <- NULL
                alternative <- 0
            }
#            anovaModel <- anovaModel0$"anovaFit"  # not used???
            boxcoxci <- c(NA, NA)
            lambda <- NA
        }

        if (type == "continuous")
        {
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
                doseFactor <- factor(cutree(hclust(dist(dose), method="average"), lenData/3))  # constructing groups containing roughly 3 observations 


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
#            if (!is.null(fct$"anovaTest")) 
            if (varPower) {testList <- mdrcVp(anovaYes = TRUE)} else {testList <- mdrcLs(anovaYes = TRUE)}
            if (!is.null(vvar)) {testList <- mdrcHetVar(anovaYes = TRUE)}
            gofTest <- testList$"gofTest"            

            if ( (!is.null(fct$"anovaYes"$"cont")) && (!is.null(testList$"anovaTest")) )
            {
                anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)
#                anovaModel0 <- (fct$"anovaTest")(anovaFormula, dset)
            } else {
                anovaModel0 <- NULL
                alternative <- 0
            }
            anovaModel <- anovaModel0$"anovaFit"


            ## Applying the Box-Cox transformation (lambda is defined here!)
            lambda <- 0

            isNumeric <- is.numeric(boxcox)
            if ( (isNumeric) || (is.logical(boxcox) && boxcox)  ) 
            {
                if (!isNumeric)
                {
                    profLik <- boxcox(anovaFormula, lambda = seq(-2.6, 2.6, 1/10), plotit=FALSE, data=dset)  # boxcox in MASS
                    maxIndex <- which.max(profLik$y)
                    lambda <- (profLik$x)[maxIndex]
                    boxcoxci <- mdrcBoxcoxCI(profLik)
                }
                if (isNumeric)
                {
                    lambda <- boxcox
                    boxcoxci <- c(NA, NA)                
                }
            } else {
                lambda <- NA
                boxcoxci <- c(NA, NA)
            }
        } 


        ## Using self starter
#        if (!is.null(fctList))
#        {
#            doseresp <- data.frame(dose, origResp)
#            
#            noSSfct <- TRUE
#            startVal <- as.vector(unlist(svList))
#            startMat <- NULL
#        }
        
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
                    startMat[i, ] <- ssfct( doseresp[(indexT1 & isfi), ] )
                } else {
                    startMat[i, ] <- rep(NA, numNames)
                }

                ## Identifying a dose response curve only consisting of control measurements
                if (sum(!is.na(startMat[i,]))==1) {upperPos <- (1:numNames)[!is.na(startMat[i,])]}
            }
            colMat <- matrix(0, numNames, numAss)
            maxParm <- rep(0, numNames)  # storing the max number of parameters
 #           startVec2 <- list()
        }
    ## Handling multi-dimensional dose data    
    } else {  
        alternative <- NULL
        anovaModel0 <- NULL
        anovaModel <- NULL
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
                indexT1 <- (assayNo==i)
                if (any(indexT1)) {startMat[i,] <- ssfct(doseresp[indexT1,])} else {startMat[i,] <- rep(NA, numNames)}

                ## Identifying a dose response curve only consisting of control measurements
                if (sum(!is.na(startMat[i,]))==1) {upperPos <- (1:numNames)[!is.na(startMat[i,])]}
            }
            colMat <- matrix(0, numNames, numAss)
            maxParm <- rep(0, numNames)  # storing the max number of parameters
#            startVec2 <- list()
        }
    }


    ## Finding parameters for the control measurements which will not be estimated
    collapseList2 <- list()
    for (i in 1:numNames)
    { 
        nccl <- ncol(collapseList[[i]])
        colNames <- colnames(collapseList[[i]])

        if ( (!is.null(cm)) && (!is.null(colNames)) ) 
        {
            accm <- as.character(cm)
            pos <- grep(accm, colNames)
            if (length(pos) == 0) 
            {
                candCol <- collapseList[[i]][, 1]
                if ( !(length(assayNoOld[candCol==1])==0) && (all(assayNoOld[candCol==1] == accm)) )
                {
                    pos <- 1  # the control measurements correspond to the "Intercept" term
                }
            }  
        } else {pos <- numeric(0)}

        if ((length(pos)>0) & !(upperPos==i)) 
        {
            collapseList2[[i]] <- as.matrix(collapseList[[i]][,-pos])  # column is removed
        } else {
            collapseList2[[i]] <- as.matrix(collapseList[[i]])  # column is kept
        } 
    }  
   
    
    ## Constructing vector of initial parameter values
    ## Filling in 0's as initial values for parameter where initial values from the self starter do not suffice
    ncclVec <- rep(0, numNames)
    startVecList <- list()
    for (i in 1:numNames)
    {
        ncclVec[i] <- ncol(collapseList2[[i]])  # ncol(as.matrix(collapseList2[[i]]))
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])
    # ncclVec and parmPos are used in 'parm2mat' function defined below


    if(!noSSfct)
    {
        nrsm <- nrow(startMat)
        for (i in 1:numNames)
        {
            sv <- rep(0, max(nrsm, ncclVec[i]))
#            indVec <- 1:min(nrsm, ncclVec[i])
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]
            
            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)
            
            startVecList[[i]] <- sv[indVec]
#            startVecList[[i]] <- sv  # indVec not used!
        }
        startVec <- unlist(startVecList)        
    } else {
        startVec <- startVal  # no checking if no self starter function is provided!!!
    }
#    print(startVec)


    ## Checking the number of start values provided
    if (!selfStart && !noSSfct) 
    {
        lenReq <- length(startVec)
        if (length(startVal) == lenReq) 
        {
            startVec <- startVal
        } else {
            stop(paste("Wrong number of initial parameter values. ", lenReq, " values should be supplied", sep=""))
        }
    }


    ## Defining function which converts parameter vector to parameter matrix            
    parm2mat <- function(parm)
    {
        parmMatrix <- matrix(0, lenData, numNames)
        for (i in 1:numNames)
        {
            parmMatrix[, i] <- collapseList2[[i]] %*% parm[parmPos[i] + 1:ncclVec[i]]
        }
        return(parmMatrix)
    }


    ## Adapting to full design matrix approach
#    if (design)
#    {
        ## Modifying initial estimates for full design matrix approach
#        print(startVec)
        
#        startMat2 <- startMat
#        startMat2[is.na(startMat2)] <- 0
#        print(startMat2)

#        mmat <- (model.matrix(~factor(assayNo) - 1))
#        print(dim(mmat))
#        print(dim(startMat2[, 1, drop=FALSE]))
#        print(as.vector(mmat%*%startMat2[, 1, drop=FALSE]))
#        print(as.vector(mmat%*%startMat2[, 2, drop=FALSE]))
#        print(as.vector(mmat%*%startMat2[, 3, drop=FALSE]))
#        print(as.vector(mmat%*%startMat2[, 4, drop=FALSE]))                
        
#        pm <- list()
#        for (i in 1:numNames)
#        {
#            clElt <- collapseList[[i]]
#            pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))%*%parm2mat2(startVec)[,i, drop=FALSE]  # ginv in MASS
#            pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))%*%parm2mat(startVec)[,i, drop=FALSE]  # ginv in MASS
#            pm[[i]] <- pm[[i]][1:ncol(collapseList[[i]])]  # maxParm[i]]
            
#            print(dim(ginv(t(clElt)%*%clElt)%*%t(clElt)))
#            indVec0 <- dim(ginv(t(clElt)%*%clElt)%*%t(clElt))[1]
#            indVec <- !is.na(startMat2[, i, drop=FALSE])
#            print(sum(indVec))
            
#            indVal <- min(c(sum(indVec), dim(clElt)[2]))
#            print(c(sum(indVec), dim(clElt)[2]))
#            print( (ginv(t(clElt)%*%clElt)%*%t(clElt))[1:indVal, ,drop=FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop=FALSE] )
#            print( ginv(t(clElt)%*%clElt)%*%t(clElt)%*%mmat%*%startMat2[, 2, drop=FALSE] )
#            print( ginv(t(clElt)%*%clElt)%*%t(clElt)%*%mmat%*%startMat2[, 3, drop=FALSE] )
#            print( ginv(t(clElt)%*%clElt)%*%t(clElt)%*%mmat%*%startMat2[, 4, drop=FALSE] )                        

            
            
#            if (ncclVec[i]<(numAss-1))  # subtracting the control
#            {
#                pm[[i]] <- startVec[parmPos[i] + 1:ncclVec[i]]
#            } else {
#            
#                mmat <- (model.matrix(~factor(assayNo) - 1))[,1:ncclVec[i], drop=FALSE]
#                curveparm <- mmat%*%(startVec[parmPos[i] + 1:ncclVec[i]])[, drop=FALSE]
#                pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))%*%curveparm
#                pm[[i]] <- pm[[i]][1:ncol(collapseList[[i]])]
#            }
            
#            if ((length(posVec[[i]])>0) && !(upperPos==i)) {pm[[i]] <- pm[[i]][-pos]}
#        } 

##        print(posVec)
##        print(upperPos)
        
#     print(parm2mat(startVec)[,3])        
#        startVec <- unlist(pm)
#        startVec <- startVec[!is.na(startVec)]

#     print(startVec)
     
     if (selfStart)
     {
#         startVec <- mdrcConvertParm(startVec, startMat, assayNo, collapseList) 
         startVec <- mdrcConvertParm(startVec, startMat, assayNo, collapseList2) 
     }
#     print(startVec)


#    ## Checking number of parameters: fct vs collapse
#    numParm2 <- 0
#    for (i in 1:numNames)
#    {
#        numParm2 <- numParm2 + ncol(as.matrix(collapseList2[[i]]))
#    }
##    print(numParm2)
#    lenReq <- length(startVec)
#    if (lenReq!=numParm2) 
#    {
#        stop(paste("Model function requires ", lenReq, " start values,\n", "        whereas 'collapse' argument requires ", numParm2, sep=""))
#    }
 

#    if (design)
#    {

#        ## Retrieving names for parameters
#        parmVecList <- list()
#        for (i in 1:numNames)
#        {
#            colNames1 <- colnames(collapseList2[[i]])
#            if (is.null(colNames1)) {parmVecList[[i]] <- paste(parNames[i],"(Intercept)",sep=":")} else { 
#                parmVecList[[i]] <- paste(parNames[i], colNames1, sep=":")}
#                
#            parmVecList[[i]] <- (parmVecList[[i]])[1:ncol(collapseList2[[i]])]  # min(maxParm[i], length(colNames1))]
#        }
#        parmVec <- unlist(parmVecList)
#        
#        parmVec2 <- parmVec
#        for (i in 1:length(parmVec))
#        {
#            pos <- regexpr("factor(collapse[, i])", parmVec[i], fixed=TRUE)
#            if (pos>0) 
#            {
#                parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos+21), sep="")
#            }
#            
#            pos <- regexpr("factor(assayNo)", parmVec[i], fixed=TRUE)
#            if (pos>0) 
#            {
#                parmVec2[i] <- paste(substring(parmVec[i],1,pos-1), substring(parmVec[i], pos+15), sep="")
#            }
#
#        }
#        parmVec <- parmVec2 

        pnList <- mdrcParNames(numNames, parNames, collapseList2)
        parmVec <- pnList[[1]]
        parmVecA <- pnList[[2]]
        parmVecB <- pnList[[3]]
        
#        print(parmVec2)
#        print(sub("factor(collapse[, i])", "", parmVec[2]))
#    }
#    }


    ## Defining non-linear function
#    indexTF <- (dose < zeroTol)  # for handling separate control measurements


#    if ((scale) && (!is.null(fct$"scaling"))) 
#    {
#        scaleList <- (fct$"scaling")(dose)
#        scaleInd <- scaleList$"scaleInd"
#
#        if (!is.na(scaleInd))
#        {        
#            scaleFct <- scaleList$"scaleFct"
#        
#            indVec <- parmPos[scaleInd] + 1:ncclVec[scaleInd]
#            startVec[indVec] <- scaleFct(startVec[indVec])
#        } else {scale <- FALSE; scaleFct <- function(x) {x}}
#    } else {scale <- FALSE; scaleFct <- function(x) {x}}
#    print(startVec)


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
            svList[[i]] <- fctList[[i]]$"ssfct"( doseresp[(indexT1 & isfi), ] )
            matList[[i]] <- c( sum(indexT1), length(svList[[i]]) )
            
            ivList2[[i]] <- match(fctList[[i]]$names, fct$names)             
        }
#        print(svList)
    
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
#            print(retVec)
            return(retVec)
        }
        
        startVec <- as.vector(unlist(svList))
    } else {
        drcFct1 <- function(dose, parm)
        {
            drcFct(dose, parm2mat(parm))
        }
    }
#print(startVec)


    if (is.null(cm)) 
    {
        multCurves <- function(dose, parm)
        {          
#           drcFct(dose, parm2mat(parm))
           drcFct1(dose, parm)  # fctList
           
#           drcFct(scaleFct(dose), parm2mat(parm))
#        
#            parmVal <- parm2mat(parm)
#
#            fctEval <- rep(0, lenData)
#            fctEval <- drcFct(scaleFct(dose), parmVal)
#
#            fctEval  # [iv]
        }
    } else {
        indexTF <- assayNoOld==cm
        indexZero <- (1:lenData)[indexTF]
        indexNZero <- (1:lenData)[!indexTF]

        drcSign <- mdrcSign(dose, resp, assayNo)
        conList <- (fct$"confct")(drcSign)
        conFct <- conList$"fct"

        multCurves<-function(dose, parm) 
        {
            parmVal <- parm2mat(parm)
                       
            fctEval <- rep(0, lenData)
            fctEval[indexZero] <- conFct(parmVal[indexZero, , drop = FALSE])
#            fctEval[indexZero] <- (fct$"confct")(parmVal[indexZero, , drop = FALSE])  # parmVal[indexZero, upperPos]
            fctEval[indexNZero] <- drcFct(dose[indexNZero], parmVal[indexNZero, , drop = FALSE])

#            assign("parmVal", parmVal, env = .GlobalEnv)  # for use with derivatives, to save calculating the same values again
#            assign("fctEval", fctEval, env = .GlobalEnv)  # for use with derivatives

            fctEval  # [iv]
        }
    }



#    if (length(indexZero)==0)  # yielding a small gain in computation time
#    {
#        multCurves<-function(dose, parm)
#        {
#            parmVal <- parm2mat(parm)
#
#            fctEval <- rep(0, lenData)
#            fctEval <- drcFct(dose, parmVal)
#
##            assign("parmVal", parmVal, env=.GlobalEnv)  # for use with derivatives, to save calculating the same values again
##            assign("fctEval", fctEval, env=.GlobalEnv)  # for use with derivatives
#
#            fctEval  # [iv]
#        }
#    } else {
#
#        multCurves<-function(dose, parm) 
#        {
#            parmVal <- parm2mat(parm)
##            print(parm2mat(parm))
##            print(parm)           
#                       
#            fctEval <- rep(0, lenData)
#            fctEval[indexZero] <- (fct$"confct")(parmVal[indexZero, ])  # parmVal[indexZero, upperPos]  # only a real assignment if 'indexTF' contains at least one TRUE
#            fctEval[indexNZero] <- drcFct(dose[indexNZero], parmVal[indexNZero, ])
#
##            assign("parmVal", parmVal, env=.GlobalEnv)  # for use with derivatives, to save calculating the same values again
##            assign("fctEval", fctEval, env=.GlobalEnv)  # for use with derivatives
#
##            print(sum((resp-fctEval)^2))
##            print(c(parm, sum(fctEval[13:15])))
#            fctEval  # [iv]
#        }
#    }

#    print(startVec)
#    print(multCurves(dose, startVec))
#    stop()


    ## Defining objective function
    robustFct <- mdrcRobust(type, robust, match.call(), lenData, length(startVec))
    

    ## Box-Cox transformation is applied
    if ( (type == "continuous") && boxcox)
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

        resp <- bcfct(resp, lambda, bcTol)

        multCurves2 <- function(dose, parm)
        {
            bcfct(multCurves(dose, parm), lambda, bcTol)
        } 
    } else {multCurves2 <- multCurves}
    

    ## Defining final objective function
    if (type == "binomial")
    {
#        iv <- !is.nan(multCurves(dose, startVec))
#        resp2 <- resp[iv]
#        weights2 <- weights[iv]
#                
#        opfct <- function(c)  # dose, resp and weights are fixed
#        {                      
#            prob <- (multCurves(dose, c))[iv]
#
#            return( -sum((resp2*weights2)*log(prob/(1-prob))+weights2*log(1-prob)) )
#        }
        estMethod <- mdrcBinomial(dose, resp, multCurves2, startVec, robustFct, weights, rmNA)    
        opfct <- estMethod$opfct 
        
        
        ## Re-fitting the ANOVA model a second time for binomial data (only if necessary)
        ## leaving out dose values 0 and infinity
        if ( (!is.null(fct$"anovaYes"$"bin")) && (!is.null(estMethod$"anovaTest2")) )
        {
            anovaModel0 <- (estMethod$"anovaTest2")(anovaFormula, dset)            
            anovaModel <- anovaModel0$"anovaFit"            
        }                            
    }
    
    if (type == "continuous")
    {
        if (varPower)
        {        
            estMethod <- mdrcVp(dose, resp, multCurves2)
            lenStartVec <- length(startVec)
            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
            parmVec <- c(parmVec, "Sigma", "Power")

        } else {
            estMethod <- mdrcLs(dose, resp, multCurves2, startVec, robustFct, weights, rmNA)
        }              
        if (!is.null(vvar))
        {
            estMethod <- mdrcHetVar(dose, resp, multCurves2, vvar)
            lenStartVec <- length(startVec)
            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
            parmVec <- c(parmVec, as.character(unique(vvar)))        
        }       
        opfct <- estMethod$opfct            


        ## Re-fitting the ANOVA model to incorporate Box-Cox transformation (if necessary)
        if (!is.na(lambda))
        {
            dset <- data.frame(dose, doseFactor, resp, assayNo, bcc)  # dataset with new resp values        
            anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)            
            anovaModel <- anovaModel0$"anovaFit"
        }
    }


    ## Re-fitting the ANOVA model to incorporate Box-Cox transformation (if necessary)
#    if (!is.na(lambda))
#    {
##        if (numAss>1) 
##        {
##            anovaModel <- lm(resp ~ offset(bcc) + factor(dose)*factor(assayNo))
##            alternative <- 2
##        } else {
##            anovaModel <- lm(resp ~ offset(bcc) + factor(dose))
##            alternative <- 1
##        }
#        dset <- data.frame(dose, doseFactor, resp, assayNo, bcc)  # dataset with new resp values        
#        anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)            
#        anovaModel <- anovaModel0$"anovaFit"
##        print(deviance(anovaModel))     
#    }
#    
#    
#    ## Re-fitting the ANOVA model a second time for binomial data (only if necessary)
#    ## leaving out dose values 0 and infinity
#    if ( (!is.null(fct$"anovaYes"$"bin")) && (!is.null(estMethod$"anovaTest2")) )
#    {
#        anovaModel0 <- (estMethod$"anovaTest2")(anovaFormula, dset)            
#        anovaModel <- anovaModel0$"anovaFit"            
#    }


    ## Defining lower and upper limits of parameters
    if (constrained)
    {
#        numParm <- table(parmVecA)
#        numParm <- numParm[match(unique(parmVecA), names(numParm))]
#        print(numParm)
#        print(ncclVec)

        if (!is.null(fct$lowerc)) 
        {
            if (!is.numeric(fct$lowerc) | !((length(fct$lowerc)==sum(ncclVec)) | (length(fct$lowerc)==numNames)))
            {
                stop("Not correct 'lowerc' argument")
            } else {
                if (length(fct$lowerc)==numNames) {lowerLimits = rep(fct$lowerc, ncclVec)} else {lowerLimits = fct$lowerc}
            } 
        }
#        print(lowerLimits) 

        if (!is.null(fct$upperc)) 
        {
            if (!is.numeric(fct$upperc) | !((length(fct$upperc)==sum(ncclVec)) | (length(fct$upperc)==numNames)))
            {
                stop("Not correct 'upperc' argument")
            } else {
                if (length(fct$upperc)==numNames) {upperLimits = rep(fct$upperc, ncclVec)} else {upperLimits = fct$upperc}
            } 
        }
#        print(upperLimits)
        
        if (all(!is.finite(lowerLimits)) & all(!is.finite(upperLimits))) {stop("No constraints are imposed via 'lowerc' and 'upperc' arguments")}
    }


    ## Optimising
##    options(warn=warnVal)
#    {if (derFlag)
#    {
#        if (constrained)
#        {
#            nlsObj <- try(optim(startVec, opfct, opfctDer, method="L-BFGS-B", lower=lowerLimits, upper=upperLimits, control=list(maxit=maxIt)), silent=TRUE)
#        } else {
#            nlsObj <- try(optim(startVec, opfct, opfctDer, method=optMethod, control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
#        }
#        options(warn=0)
#        
#        if (!inherits(nlsObj, "try-error")) 
#        {
#            nlsFit <- nlsObj
#        } else {
##            stop("Convergence failed")
#            warning("Convergence failed. The model was not fitted!", call.=FALSE)
#
#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}
#
#            return(list(call=callDetail, parNames=parmVec, startVal=startVec))
#        }
#        nlsFit$hessian <- opfctDer2(parmVal)
#
#    } else {  # in case no derivatives are used
#
#        if (constrained)
#        {
#        nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method="L-BFGS-B", lower=lowerLimits, upper=upperLimits, control=list(maxit=maxIt)), silent=TRUE)
#        } else {
##            print(startVec)
#            psVec <- abs(startVec)
#            psVec[psVec<1e-4] <- 1
##            psVec <- rep(1, length(startVec))
#            nlsObj <- try(optim(startVec, opfct, hessian=TRUE, method=optMethod, control=list(maxit=maxIt, reltol=relTol, parscale=psVec)), silent=TRUE)
#        }
#        options(warn=0)
#        
#        if (!inherits(nlsObj, "try-error")) 
#        {
#            nlsFit <- nlsObj
#        } else {
#            if (errorMessage) {stop("Convergence failed")} else {warning("Convergence failed. The model was not fitted!", call.=FALSE)}
#
#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}
#
#            return(list(call=callDetail, parNames=parmVec, startVal=startVec))
#
#        }
#    }}

##    options(warn=0)
##    return(nlsFit)


    ## Manipulating before optimisation
    scaleInd <- fct$"scaleInd"
    if (!is.null(scaleInd))
    {
        if (!is.null(fctList))
        {
            parmInd <- rep(0, numAss)
            for (i in 1:numAss)
            {
                parmInd[i] <- fctList[[i]]$"scaleInd"
            }
            parmInd <- cumsum(parmInd)
        } else {
            parmInd <- parmPos[scaleInd] + 1:ncclVec[scaleInd]        
        }
    
        scaleFct <- mdrcScaleDose(dose, maxDose)        
        dose <- scaleFct(dose)
        startVec[parmInd] <- scaleFct(startVec[parmInd])
    }
#print(startVec)
#print(multCurves2(dose, startVec))

    ## Optimising
    nlsFit <- mdrcOpt(opfct, startVec, optMethod, derFlag, constrained, warnVal, upperLimits, lowerLimits, errorMessage, maxIt, relTol) 
    nlsFit$ofvalue <- nlsFit$value


    ## Manipulating after optimisation
    if (!is.null(scaleInd))
    {
        dose <- scaleFct(dose, down = FALSE)
#        parmInd <- parmPos[scaleInd] + 1:ncclVec[scaleInd]
        startVec[parmInd] <- scaleFct(startVec[parmInd], down = FALSE)
        nlsFit$par[parmInd] <- scaleFct(nlsFit$par[parmInd], down = FALSE)
#        nlsFit$hessian[, parmInd] <- scaleFct(nlsFit$hessian[, parmInd])
#        nlsFit$hessian[parmInd, ] <- scaleFct(nlsFit$hessian[parmInd, ])

        scaleFct2 <- function(hessian) 
                     {
                         newHessian <- hessian
                         newHessian[, parmInd] <- scaleFct(newHessian[, parmInd], down = FALSE)
                         newHessian[parmInd, ] <- scaleFct(newHessian[parmInd, ], down = FALSE)
                         return(newHessian)
                     }                
    } else {scaleFct2 <- function(x) {x}}
#    print(dose)
#    print(startVec)


#    ## Scaling estimates back
#    if (scale)
#    {
#        nlsFit$par[indVec] <- scaleFct(nlsFit$par[indVec], down=FALSE)
#        
##        nlsFit$hessian[indVec]
#    }


    ## Handling variance parameters
#    vpType <- NULL    
#    vpIndex <- NULL
    varParm <- NULL
        
    if (varPower)
    {
#        vpType <- "varPower"
#        vpIndex <- 1:lenStartVec  # (lenStartVec+1):length(nlsFit$par)
        varParm <- list(type = "varPower", index = 1:lenStartVec)        
    }
# else {
#        varParmType <- NULL    
#        varParmIndex <- NULL
#    }
    if (!is.null(vvar))
    {
#        vpType <- "hetvar"
#        vpIndex <- 1:lenStartVec  # (lenStartVec+1):length(nlsFit$par)    
        varParm <- list(type = "hetvar", index = 1:lenStartVec)
    }
#    varParm <- list(type=vpType, index=vpIndex)


    # Testing against the ANOVA (F-test)
#    if (!is.null(anovaModel) && anovaFlag)
#    {
#        anovaSS <- deviance(anovaModel)
#        anovaDF <- df.residual(anovaModel)
        nlsSS <- nlsFit$value
        nlsDF <- lenData - length(startVec)
#        ftest <- (nlsSS-anovaSS)/(nlsDF-anovaDF)/(anovaSS/anovaDF)
#        pVal2 <- pf(ftest,nlsDF-anovaDF,anovaDF,lower.tail=FALSE)
#    } else {
#        anovaSS <- NA
#        anovaDF <- NA
#        nlsSS <- NA
#        nlsDF <- NA
#        ftest <- NA
#        pVal2 <- NA
#    }


    ## Constructing a plot function
        
    ## Picking parameter estimates for each curve. Does only work for factors not changing within a curve!
#    if (exists("assayNZero")) {iVec <- (1:numAss)[-assayNZero]} else {iVec <- 1:numAss}
    if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames==cm)]} else {iVec <- 1:numAss}
    
    pickCurve <- rep(0, length(iVec))
    for (i in iVec)
    {
       pickCurve[i] <- (1:lenData)[assayNo==i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)
#    parmMat[iVec,] <- parm2mat(nlsFit$par)[pickCurve,]

    fixedParm <- (estMethod$"parmfct")(nlsFit)
    parmMat[iVec,] <- parm2mat(fixedParm)[pickCurve,]

    if(!is.null(fctList))
    {
         parmMat <- matrix(NA, numAss, numNames)
         for (i in 1:numAss)
         {
#             print(ivList2[[i]])
#             print(fixedParm[(posVec[i]+1):posVec[i+1]])
             parmMat[i, ivList2[[i]]] <- fixedParm[(posVec[i]+1):posVec[i+1]]
         }
    }    
    

#    if (exists("assayNZero")) 
#    {
#        parmMat[-iVec, upperPos] <- parm2mat(nlsFit$par)[assayNo==assayNZero, ][1, upperPos]  # all rows are identical ... taking the first one
#        parmMat[-iVec, upperPos] <- parm2mat(fixedParm)[assayNoOld==assayNZero, ][1, upperPos]
#    }
    
    if (!is.null(cm))
    {
        conPos <- conList$"pos"
#        print(conPos)
        parmMat[-iVec, conPos] <- parm2mat(fixedParm)[assayNoOld==cm, , drop=FALSE][1, conPos]  # 1: simply picking the first row
    }
    rownames(parmMat) <- assayNames
#    print(parmMat)


    pmFct <- function(fixedParm)
    {
        if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames == cm)]} else {iVec <- 1:numAss}
    
#        pickCurve <- rep(0, length(iVec))
#        for (i in iVec)
#        {
#            pickCurve[i] <- (1:lenData)[assayNo==i][1]
#        }
#        parmMat <- matrix(NA, numAss, numNames)
#
##        fixedParm <- (estMethod$"parmfct")(nlsFit)
#        parmMat[iVec,] <- parm2mat(fixedParm)[pickCurve,]
    
        if (!is.null(cm))
        {
            conPos <- conList$"pos"
            parmMat[-iVec, conPos] <- parm2mat(fixedParm)[assayNoOld==cm, , drop=FALSE][1, conPos]  # 1: simply picking the first row
        }
        rownames(parmMat) <- assayNames
    
        return(parmMat)
    }
    parmMat <- pmFct( (estMethod$"parmfct")(nlsFit) )


    ## Constructing design matrix allowing calculations for each curve
    colPos <- 1
    rowPos <- 1
#    Xmat <- matrix(0, numAss*numNames, length(nlsFit$par))
    Xmat <- matrix(0, numAss*numNames, length(fixedParm))


    if (!is.null(fctList)) {omitList <- list()}
    for (i in 1:numNames)
    {
        indVec <- iVec
        lenIV <- length(indVec)

        nccl <- ncol(collapseList2[[i]])  # min(maxParm[i], ncol(collapseList2[[i]]))        

        XmatPart <- matrix(0, lenIV, nccl)
        k <- 1
        if (!is.null(fctList)) {omitVec <- rep(TRUE, lenIV)}
        for (j in indVec)
        {
            if (!is.null(fctList))
            {
                parPresent <- !is.na(match(i, ivList2[[j]]))
                omitVec[k] <- parPresent
            }
            
#            if (!is.na(parPresent))
#            {
                XmatPart[k, ] <- (collapseList2[[i]])[(1:lenData)[assayNo == j][1], 1:nccl]
#            }  # else {
#                print(j)
#                XmatPart <- XmatPart[-j, , drop = FALSE]
#            }
            k <- k + 1
 #           print(XmatPart)
        }
        if (!is.null(fctList))
        {
            XmatPart <- XmatPart[omitVec, , drop = FALSE]
            nccl <- nccl - sum(!omitVec)
            omitList[[i]] <- omitVec
        }


        Xmat[rowPos:(rowPos+lenIV-1), colPos:(colPos+nccl-1)] <- XmatPart
        colPos <- colPos + nccl
        rowPos <- rowPos + lenIV
    }
    Xmat <- Xmat[1:(rowPos-1), 1:(colPos-1)]


    ## Defining the plot function
#    plotFct <- function(dose)
#    {
#        if (doseDim==1) {lenPts<-length(dose)} else {lenPts<-nrow(dose)}
#
#        curvePts<-matrix(NA, lenPts, numAss)
#        for (i in 1:numAss)
#        {
#            if (!is.null(fctList)) 
#            {
#                drcFct <- fctList[[i]]$"fct"
#                numNames <- matList[[i]][2]    
#            }
#            print(numNames)
#        
#            if (i%in%iVec)
#            {
##                parmChosen <- parmMat[i, ]
#                parmChosen <- parmMat[i, complete.cases(parmMat[i, ])]  # removing NAs
#                
#                
#                parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow=TRUE)
#                #if (exists("assayNZero") && assayNZero==i) {curvePts[,i] <- NA} else {curvePts[,i] <- drcFct(dose, parmMat2)}
#                curvePts[,i] <- drcFct(dose, parmMat2)
#            } else { curvePts[,i] <- rep(NA, lenPts)}
#        }
#        return(curvePts)
#    }

    
    pfFct <- function(parmMat)
    {
        plotFct <- function(dose)
        {
            if (doseDim == 1) {lenPts <- length(dose)} else {lenPts <- nrow(dose)}

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
                    
                    parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow=TRUE)
                    curvePts[,i] <- drcFct(dose, parmMat2)
                } else { curvePts[,i] <- rep(NA, lenPts)}
            }
            return(curvePts)
        }
    
        return(plotFct)    
    }
    plotFct <- pfFct(parmMat)


    ## Computation of fitted values and residuals
#    if (boxcox) 
#    {           
##        predVec <- bcfct(multCurves(dose, nlsFit$par), lambda, bcTol)    
#        predVec <- bcfct(multCurves(dose, fixedParm), lambda, bcTol)    
#    } else {
##        predVec <- multCurves(dose, nlsFit$par)
#        predVec <- multCurves(dose, fixedParm)
#    }

#print(fixedParm)
#print(dose)

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
    if (type=="binomial")
    {
        nlsFit$value <- nlsDF/2
    }

    
    ## Adding meaningful names for robust methods
    robust <- switch(robust, median="median", trimmed="metric trimming", tukey="Tukey's biweight", 
                             winsor="metric Winsorizing", lms="least median of squares",
                             lts="least trimmed squares")


    ## Creating a kind of design matrix
#    collapse <- as.matrix(collapse)
#    colnames(collapse) <- parNames


    ## The summary vector
    sumVec <- c(lambda, NA, NA, NA, nlsSS, nlsDF, lenData, alternative)


    sumList <- list(lenData = lenData, alternative = alternative, df.residual = lenData - length(startVec))


    ## The function call
    callDetail <- match.call()
    if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}

    ## The data set
    if (!is.null(logDose)) {dose <- origDose}
    dataSet <- data.frame(dose, origResp, assayNo, assayNoOld, weights)
    names(dataSet) <- c(varNames, anName, anName, "weights")


    ## Detaching packages used 
#    detach(package:MASS)


    ## Box-Cox information
    bcVec <- c(lambda, boxcoxci)
    if (all(is.na(bcVec))) {bcVec <- NULL}
    if (!is.null(bcVec)) {bcVec <- c(bcVec, bcAdd)}


    ## Evaluating goodness-of-fit test
    if (!is.null(gofTest)) {gofTest <- gofTest(resp, weights, predVec, sumList$"df.residual")}


    ## Defining function that returns indices for parameters and var-cov matrix for a particular curve
#    getCurveInd <- function(curve)
#    {
#
#        if (!curve%in%unique(assayNoOld)) {stop(paste("Specified curve", curve, "not in dataset"))}
#
#        firstInd <- (tapply(1:lenData, assayNo, function(x){x[1]}))[unique(assayNoOld)%in%curve]
#
#        indVec <- list()
#        parVec <- rep(0, numNames)
#        for (i in 1:numNames)
#        {
#            indVec[[i]] <- collapseList[[i]][firstInd,]
#            parVec[i] <- ncol(collapseList[[i]])
#        }
#        parVec <- cumsum(parVec)
#        retVec <- as.vector(unlist(indVec))
#
#        lenRV <- length(retVec)
#        transMat <- matrix(0, numNames, lenRV)
#
#        for (i in 1:numNames)
#        {
#            indVec2 <- max(parVec[i-1]+1,1):parVec[i]
#            transMat[i, (1:lenRV)[indVec2]] <- retVec[indVec2]
#        }
#        return(list(retVec, transMat))
#    }



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
 

    ## Returning the fit
    returnList <- list(varParm, nlsFit, list(plotFct, logDose), sumVec, startVec, list(parmVec, parmVecA, parmVecB), 
    diagMat, callDetail, dataSet, t(parmMat), fct, Xmat, robust, type, bcVec, estMethod, lenData-length(startVec), 
    anovaModel0, gofTest, sumList, scaleFct2, pmFct, pfFct)
    names(returnList) <- c("varParm", "fit", "curve", "summary", "startVal", "parNames", "predres", "call", "data", 
    "parmMat", "fct", "transformation", "robust", "type", "boxcox", "estMethod", "df.residual", "anova", "gofTest", 
    "sumList", "scaleFct", "pmFct", "pfFct")
    class(returnList) <- c("drc", class(fct))

    return(returnList)
}
