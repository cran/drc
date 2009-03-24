\name{vcov.drc}

\alias{vcov.drc}

\title{Calculating variance-covariance matrix for objects of class 'drc'}

\description{
  'vcov' returns the estimated variance-covariance matrix for the parameters in the non-linear function.
}

\usage{
  \method{vcov}{drc}(object, ..., corr = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{corr}{logical. If TRUE a correlation matrix is returned.}
}

\value{
  A matrix of estimated variances and covariances.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
vcov(ryegrass.m1)
vcov(ryegrass.m1, corr = TRUE)

}
\keyword{models}
\keyword{nonlinear}
