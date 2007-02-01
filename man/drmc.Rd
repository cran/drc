\name{drmc}

\alias{drmc}

\title{Sets control arguments}

\description{
  Set control arguments in the control argument in the function 'drm'.
}

\usage{
  drmc(constr = FALSE, errorm = TRUE, maxDose = 1e-1, maxIt = 500, method="BFGS", 
  noMessage = FALSE, relTol = 1e-07, rmNA=FALSE, useD = FALSE, warnVal = -1, zeroTol=0)
}

\arguments{
  \item{constr}{logical. If TRUE optimisation is constrained, only yielding non-negative parameters.}
  \item{errorm}{logical specifying whether failed convergence in \code{\link{drm}} should result 
    in an error or only a warning.}
  \item{maxDose}{numeric. The maximum dose allowed if scaling of dose is to be used.}
  \item{maxIt}{numeric. The maximum number of iterations in the optimisation procedure.}
  \item{method}{character string. The method used in the optimisation procedure. 
    See \code{\link{optim}} for available methods.}
  \item{noMessage}{logical, specifying whether or not messages should be displayed.}
  \item{relTol}{numeric. The relative tolerance in the optimisation procedure.}
  \item{rmNA}{logical. Should NAs be removed from sum of squares used for estimation? 
    Default is FALSE (not removed).}
  \item{useD}{logical. If TRUE derivatives are used for estimation (if available).}  
  \item{warnVal}{numeric. If equal to 0 then the warnings are stored and displayed at the end. 
    See under 'warn' in \code{\link{options}}. The default results in suppression of warnings.}
  \item{zeroTol}{numeric. The threshold below which doses are considered as control.}
}

\value{
  A list with 8 components, one for each of the above arguments.
}

\author{Christian Ritz}

\note{
  The use of a non-zero constant \code{bcAdd} may in some cases 
  make it more difficult to obtain convergence of the estimation procedure.}

%\seealso{}

\examples{

### Displaying the default settings
drmc()  

### Using 'method' argument
model1 <- drm(ryegrass, fct = LL.4())

model2 <- drm(ryegrass, fct = LL.4(),
control = drmc(method = "Nelder-Mead"))

}
\keyword{models}
\keyword{nonlinear}
