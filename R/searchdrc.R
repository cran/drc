"searchdrc" <- function(object, which, range, len=50)
{

   sv <- object$startVal
 
   parNames <- object$parNames
   whichInd <- regexpr(paste("^",which,":",sep=""), parNames)
   whichInd <- ((1:length(parNames))[whichInd>0])[1]
   
   if (length(whichInd)<1) {stop(paste("No such parameter ", which, sep=""))}

   found <- FALSE
   for (i in seq(range[1], range[2], length.out=len))
   {

       sv[whichInd] <- i 

       options(warn=-1)
       modeldrc <- try(update(object, startVal=sv, control=mdControl(errorm=TRUE)), silent=TRUE)
       if (!inherits(modeldrc,"try-error")) {found <- TRUE; break}
       options(warn=0)
   }
   
   if (found) 
   {
      return(modeldrc)
   } else {
      warning("Convergence failed.", call. = FALSE)
   }
}
