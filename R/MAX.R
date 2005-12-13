"MAX" <-
function(obj, upper=NULL)
{
#    typeStr <- "ED"
#
#
#    ## Finding super class
#    EDinfo <- partIn(class(obj), typeStr)
#
#    EDstr <- EDinfo[[1]]
#    EDfct <- EDinfo[[2]]
#
#    inOut <- iofct(typeStr, EDstr)
#    in1fct <- inOut[[1]]
#    in2fct <- inOut[[2]]
#    outfct <- inOut[[3]]


    ## Checking class of 'obj'
    MAXlist <- obj[[11]]$"maxfct"
    if (is.null(MAXlist)) {stop("No method available")}


    ## Retrieving relevant quantities
    assayNo<-obj[[9]][,3]
    numAss<-length(unique(assayNo))

    parmMat <- obj[[10]]
    sumObj<-summary(obj)
    resVar<-sumObj[[1]]
    varMat<-obj[[12]]%*%sumObj[[2]]%*%t(obj[[12]])
    parm<-c((sumObj[[3]])[,1])
#    strParm <- (unlist(strsplit(obj[[6]], ":")))[(1:length(obj[[6]]))*2] 
    strParm <- unique(obj[[9]][, ncol(obj[[9]]) - 1])  # second last column contains original curve levels
#    strParm <- strParm[apply(parmMat, 2, function(x){!any(is.na(x))})]

#    ncPM <- ncol(parmMat)
#    naVec <- rep(NA, ncPM)
#    for (i in 1:ncPM)
#    {
#        if (any(is.na(parmMat[,i]))) {naVec[i] <- i}
#    }


    ## Creating an index matrix
#    asVec1 <- c(t(parmMat))
#    notNA <- !is.na(asVec1)
#    asVec2 <- asVec1[notNA]
#    asVec1[notNA] <- match(asVec2, unique(asVec2))
#    indexMat <- matrix(asVec1, nrow(parmMat), ncol(parmMat), byrow=TRUE)

#    lenNV <- length(naVec[!is.na(naVec)])
#    if (lenNV>0)
#    {
#        parmMat <- parmMat[,-naVec[!is.na(naVec)]]
#        indexMat <- indexMat[,-naVec[!is.na(naVec)]]
#    }
    
    ncPM <- ncol(obj[[10]])
    nrPM <- nrow(obj[[10]])
    
    options(warn=-1)  # to avoid warnings when filling in matrix with elements in excess in the vector 
    indexMat <- t(matrix(NA, nrPM, ncPM))
    indexMat[!is.na(t(obj[[10]]))] <- 1:(nrPM*ncPM)
    indexMat <- t(indexMat)
    options(warn=0)
#    print(indexMat)
    
    naVec <- rep(FALSE, ncPM)
    for (i in 1:ncPM)
    {
        naVec[i] <- any(is.na(indexMat[,i]))
    }
    indexMat <- indexMat[, !naVec, drop=FALSE]
    parmMat <- parmMat[, !naVec, drop=FALSE] 
#    print(indexMat)   

    ## Finding out which parameter occurs most times; this determines the number of ED values
#    maxIndex <- 0
#    maxParm <- 0
#    
#    print(EDlist(parmMat[,1], 50, upper))
#    indexVec <- (1:nrow(parmMat))[(EDlist(parmMat[,1], 50, upper))[[2]]<1e-10]
#    print(indexVec)
#    for (i in indexVec)
#    {
#        PM <- parmMat[i,]
#        lenPM <- length(unique(PM))
#        if (lenPM > maxParm) {maxIndex <- match(unique(PM),PM); maxParm <- i}
#    }
#
#    nCol <- ncol(parmMat)
#    indexVec <- 1:nCol
#    for (i in 1:nCol)
#    {
#        if (any(is.na(parmMat[,i]))) {indexVec[i] <- NA}
#    }
#    indexVec <- indexVec[!is.na(indexVec)]
#    print(indexVec)
    indexVec <- 1:ncol(indexMat)


    ## Calculating ED values    
    lenEB <- length(indexVec)    
    dimNames<-rep("", lenEB)
    MAXmat<-matrix(0, lenEB, 2)

    
#    rowIndex <- 1
#    for (i in maxIndex)
    for (i in indexVec)
    {
        parmInd <- indexMat[,i]
        varCov <- varMat[parmInd, parmInd]
        parmChosen <- parmMat[,i]

#        print(parmChosen)
        MAXmat[i, ] <- MAXlist(parmChosen, upper)
        dimNames[i] <- strParm[i]
    }

    dimnames(MAXmat)<-list(dimNames,c("Dose","Maximum"))
    printCoefmat(MAXmat)
    invisible(MAXmat)    
#    return(EDmat)
}
