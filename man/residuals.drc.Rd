\name{residuals.drc}

\alias{residuals.drc}

\title{Extracting residuals from model}

\description{
  'residuals' extracts residuals from an object of class 'drc'.
}

\usage{

  \method{residuals}{drc}(object, rstandard = FALSE, ...)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{rstandard}{logical. If TRUE standardised residuals are returned.}
  \item{...}{additional arguments.}
}

\value{
  Raw or standardised residuals extracted from 'object'.
}

\note{
  The argument 'rstandard' is available for least squares estimation
  with or without Box-Cox transformation or variance as a power of the
  mean.
}

\author{Christian Ritz}

\examples{

model1 <- multdrc(FA)

## Displaying the residual plot (raw residuals)
plot(fitted(model1), residuals(model1))

## Using the standardised residuals
plot(fitted(model1), residuals(model1, rstandard = TRUE))

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
