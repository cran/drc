"fitted.drc" <-
function(object, ...)
{
    if (inherits(object, "bindrc"))
    {
        olddata <- object$data
    
        newdata <- data.frame(dose=olddata$dose, int=olddata$int, slope=olddata$slope)
        
        if (!is.null(olddata$low)) {newdata <- data.frame(newdata, low = olddata$low)}
        if (!is.null(olddata$up)) {newdata <- data.frame(newdata, up = olddata$up)}
    
        return(predict(object, newdata))
    } else {

        return(object[[7]][,1])
    }
}
