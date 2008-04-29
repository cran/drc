\name{mdControl}

\alias{mdControl}

\title{Sets control arguments}

\description{
  Set control arguments in the control argument in the function 'multdrc'.
}

\usage{
mdControl(bcAdd=0, constr = FALSE, errorm = TRUE, maxDose = 1e-1, maxIt = 500, method="BFGS", 
          noMessage = FALSE, relTol = 1e-07, rmNA=FALSE, warnVal = -1, zeroTol=0)
}

\arguments{
  \item{bcAdd}{numeric. The constant to be added in Box-Cox transformation.}
  \item{constr}{logical. If TRUE optimisation is constrained, only yielding non-negative parameters.}
  \item{errorm}{logical specifying whether failed convergence in \code{\link{multdrc}} should result in an error or only a warning.}
  \item{maxDose}{numeric. The maximum dose allowed if scaling of dose is to be used.}
  \item{maxIt}{numeric. The maximum number of iterations in the optimisation procedure.}
  \item{method}{character string. The method used in the optimisation procedure. See \code{\link{optim}} for available methods.}
  \item{noMessage}{logical, specifying whether or not messages should be displayed.}
  \item{relTol}{numeric. The relative tolerance in the optimisation procedure.}
  \item{rmNA}{logical. Should NAs be removed from sum of squares used for estimation? Default is FALSE (not removed).}
  \item{warnVal}{numeric. If equal to 0 then the warnings are stored and displayed at the end. See under 'warn' in \code{\link{options}}.
                 The default results in suppression of warnings.}
  \item{zeroTol}{numeric. The threshold below which doses are considered as control.}
}

\value{
  A list with 8 components, one for each of the above arguments.
}

\author{Christian Ritz}

\note{The use of a non-zero constant \code{bcAdd} may in some cases make it more difficult to obtain convergence of the estimation procedure.}

\seealso{See the documentation on the data set \code{\link{spinach}}.}

\examples{

### Displaying the default settings
mdControl()  

### Using 'method' argument
ryegrass.model1 <- multdrc(ryegrass)

ryegrass.model2 <- multdrc(ryegrass, varPower = TRUE, startVal = coef(ryegrass.model1),
control = mdControl(method = "Nelder-Mead"))

### Using 'relTol' argument to improve accuracy

## Fitting a model to the PestSci dataset
spinach.model1 <- multdrc(SLOPE~DOSE, CURVE, data = spinach)

## Fitting the same model in a different parametrisation
spinach.model2 <- multdrc(SLOPE~DOSE, CURVE, data = spinach,
collapse = list(~factor(CURVE),~factor(CURVE),~factor(CURVE),~factor(CURVE)))

## Reducing the relative tolerance
spinach.model3 <- multdrc(SLOPE~DOSE, CURVE, data = spinach,
collapse = list(~factor(CURVE),~factor(CURVE),~factor(CURVE),~factor(CURVE)),
control = mdControl(relTol = 1e-12))

## Comparing the parameter estimates
coef(spinach.model1)
coef(spinach.model2)
coef(spinach.model3)  # model3 is closer to model1 than model2 is 

}
\keyword{models}
\keyword{nonlinear}
