\name{SI}
\alias{SI}

\title{Comparing selectivity indices across curves}

\description{
  'SI' compares selectivity indices for arbitrary dosage across curves.
  The selectivity is the ratio between effective dosages from different curves. 
}

\usage{
SI(obj, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, ...)
}

\arguments{
  \item{obj}{an object of class 'drc'.}
  \item{percVec}{a numeric vector of dosage values.}
  \item{compMatch}{an optional character vector of names of assays to be compared. If not specified all comparisons are supplied.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{reverse}{logical. If TRUE the order of comparison of two curves is reversed.}
  \item{...}{additional arguments to the function doing the calculations. 
  For instance Upper limit for bisection method. 
  Need to be larger than ED values to be calculated in ratio.}  
}

\details{  
  For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper' argument
  or the 'interval' argument. The 'upper' argument specifies the upper limit of the bisection method. 
  The upper limits needs to be larger than the EDx level to be calculated. The default limit is 1000.
  The 'interval' argument should specify a rough interval in which the dose
  yielding the maximum hormetical response lies. The default interval is 'c(0.001, 1000)'. Notice that the lower limit
  should not be set to 0 (use something like 1e-3, 1e-6, ...).
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, 
  t-statistics for testing indices equal to 1 and the corresponding p-values.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}

\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\seealso{The related function \code{\link{ED}}.}

\examples{

model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)

SI(model1, c(50,50))
SI(model1, c(10,50))
SI(model1, c(10,50), reverse=TRUE)

rm(model1)

}

\keyword{models}
\keyword{nonlinear}
