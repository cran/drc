\name{ED}
\alias{ED}

\title{Estimating effective dosage}

\description{
  'ED' estimates effective dosage values (EDx) using a single model fit.
}

\usage{
  ED(obj, percVec, od = FALSE, ...)
}

\arguments{
  \item{obj}{an object of class 'drc'.}
  \item{percVec}{a numeric vector containing the dosage values.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{...}{additional arguments for the function doing the calculations.}
}

\details{
  This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.

  For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper' argument
  or the 'interval' argument. 
  The 'upper' argument specifies the upper limit of the bisection method. 
  The upper limit needs to be larger than the EDx level to be calculated. The default limit is 1000, but
  this may need to be increased.
  The 'interval' argument should specify a rough interval in which the dose
  yielding the maximum hormetical response lies. The default interval is 'c(0.001, 1000)'. Notice that the lower limit
  should not be set to 0 (use something like 1e-3, 1e-6, ...).
}

\value{
  A matrix with two columns: one containing the estimates and one containing the corresponding estimated standard errors.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\seealso{The related function \code{\link{SI}}.}

\examples{

model1<-multdrc(FA)

ED(model1, c(10,50,90))[,1]   
## first column: the estimates of ED10, ED50 and ED90
## second column: the estimated standard errors 

rm(model1)

}

\keyword{models}
\keyword{nonlinear}
