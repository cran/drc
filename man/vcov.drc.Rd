\name{vcov.drc}

\alias{vcov.drc}

\title{Calculating Variance-Covariance Matrix for objects of class 'drc'}

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

model1 <- multdrc(ryegrass)
vcov(model1)
vcov(model1, corr = TRUE)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
