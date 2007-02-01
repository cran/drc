\name{logLik.drc}

\alias{logLik.drc}

\title{Extracting the log likelihood}

\description{
  'loglik' extracts the value of the log likelihood function evaluated at the parameter estimates.
}
\usage{
  \method{logLik}{drc}(object, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
   The evaluated log likelihood as a numeric value and the corresponding degrees of freedom as an attribute.
}

\note{
    The value of the log likelihood could be used to compare fits of the same data and model for different starting values. Larger is better.
}

\author{Christian Ritz}

\examples{

model1 <- multdrc(ryegrass) 
logLik(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
