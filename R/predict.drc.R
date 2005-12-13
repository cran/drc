"predict.drc" <- function(object, newdata, ...)
{
    if (inherits(object, "bindrc"))
    {

    dose <- newdata$dose
    if (is.null(dose)) {stop("No dose values provided")}
    lenData <- length(dose)
    
    if (object$"logTrans") {ldose <- log(dose)} else {ldose <- dose}

    int <- newdata$int
    if (is.null(int)) {intVal <- rep(object$predict$int[1], lenData)} else {intVal <- object$predict$int[as.character(int)]}
#    print(int)

    slope <- newdata$slope
    if (is.null(slope)) {intVal <- rep(object$predict$slope[1], lenData)} else {slopeVal <- object$predict$slope[as.character(slope)]}

    lower <- newdata$low
    if (is.null(lower)) {lowVal <- rep(0, lenData)} else {lowVal <- object$predict$low[as.character(lower)]}
    if ( (is.null(lowVal)) || (all(is.na(lowVal)))) {lowVal <- rep(0, lenData)}

    upper <- newdata$up
    if (is.null(upper)) {upVal <- rep(1, lenData)} else {upVal <- object$predict$up[as.character(upper)]}
    if ( (is.null(upVal)) || (all(is.na(upVal))) ) {upVal <- rep(1, lenData)}

#    intVal <- object$predict$int[int] 
#    slopeVal <- object$predict$slope[slope] 
#    print(intVal) 
#    print(slopeVal) 
#    print(lowVal)

    link <- object$link  # this is necessary!
    predVal <- lowVal + (upVal - lowVal)*binomial(link)$linkinv(intVal + slopeVal*ldose)

    return(as.vector(predVal))
    } else {stop("No 'predict' method available")}
}
