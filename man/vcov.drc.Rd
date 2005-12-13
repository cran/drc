\name{vcov.drc}
\alias{vcov.drc}
\title{Calculating Variance-Covariance Matrix for objects of class 'drc'}
\description{
  'vcov' returns the estimated variance-covariance matrix for the parameters in the non-linear function.
}
\usage{
vcov.drc(object, ...)
}
\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}
\value{
  A matrix of estimated variances and covariances.
}

\author{Christian Ritz}

\examples{

model1 <- multdrc(FA)
vcov(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
