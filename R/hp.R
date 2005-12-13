"hp" <- function(lowerc=c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), upperc=c(Inf, Inf, Inf, Inf, Inf, Inf), fixed=c(NA, NA, NA, NA, NA, NA), 
                 names=c("b","c","d","e","f","g"), useDer=FALSE)
{
    ## Checking arguments
    numParm <- 6
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
    if (useDer) {stop("Derivatives not available")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
#        loge <- -parm[, 6]*log((1/parm[, 4])^(1/parm[, 6]) + (1/parm[, 5])^(1/parm[, 6]))
        loge <- -parm[, 6]*log((parm[, 4])^(1/parm[, 6]) + (parm[, 5])^(1/parm[, 6]))
        parmMat[, 2]+(parmMat[, 3]-parmMat[, 2])/(1+exp(parmMat[, 1]*(log(dose)-loge)))
    }


    ## Defining value for control measurements (dose=0)
    confct <- function(drcSign)
    {
        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
        confct2 <- function(parm)
        { 
            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
            parmMat[, notFixed] <- parm
            parmMat[, conPos]
        }
        return(list(pos=conPos, fct=confct2))
    }


    ## Defining flag to indicate if more general ANOVA model
    anovaYes <- FALSE


    ## Defining the self starter function
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

        startVal[3] <- max(resp3)+0.001  # the d parameter
        startVal[2] <- min(resp3)-0.001  # the c parameter
        startVal[5] <- 1  # better choice may be possible!        
#        startVal[!notFixed] <- fixed[!notFixed] 
        
        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[3]-resp3)/(resp3-startVal[2]+0.001))  # 0.001 to avoid 0 in the denominator
        logitFit <- lm(logitTrans~log(dose3))
        startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
        
        startVal[5] <- startVal[4]
        startVal[6] <- 1

        return(startVal[notFixed])
    }
    ssApList <- list(mean, mean, mean, function(x) x[uniqueNames=="100"], function(x) x[uniqueNames=="0"], mean)
 
   
    ## Defining names
    names <- names[notFixed]


    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL


    ## Limits
    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- NULL


    ## Defining the SI function
    sifct <- NULL


    returnList <- list(fct=fct, confct=confct, anovaYes=anovaYes, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2, 
                       lowerc=lowerLimits, upperc=upperLimits, edfct=edfct, sifct=sifct)
                       
    class(returnList) <- "Hewlett-Plackett"
    invisible(returnList)
}
# model2<-multdrc(rgr~dose,pct,data=TM,collapse=list(~factor(pct),~1,~1,~I(pct/100),~I(1-pct/100),~1),fct=hp(),startVal=sv)
# model2<-multdrc(rgr~dose,pct,data=TM,collapse=list(~factor(pct),~1,~1,~I(pct/100)-1,~I(1-pct/100)-1,~1),fct=hp(),startVal=sv)

#> sv
#   b:(Intercept)  b:factor(pct)17  b:factor(pct)33  b:factor(pct)50 
#      1.41444870      -0.09480835       0.11094586      -0.02216643 
# b:factor(pct)67  b:factor(pct)83 b:factor(pct)100    c:(Intercept) 
#     -0.05609425      -0.63248955      -1.22583086      -0.04846963 
#   d:(Intercept)         e:I(pct)    f:(Intercept)    g:(Intercept) 
#      0.29084899      10.08429917    6621.70600293      -0.01243392

#sv <- c(1.41444870, -0.09480835, 0.11094586, -0.02216643,-0.05609425, -0.63248955, -1.22583086, -0.04846963, 0.29084899, 10.08429917, 6621.70600293,-0.5)
#model2<-multdrc(rgr~dose,pct,data=TM,collapse=list(~factor(pct),~1,~1,~I(pct/100)-1,~I(1-pct/100)-1,~1),fct=hp(),startVal=sv)
