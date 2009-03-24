\name{residuals.drc}

\alias{residuals.drc}

\title{Extracting residuals from model}

\description{
  'residuals' extracts residuals from an object of class 'drc'.
}

\usage{

  \method{residuals}{drc}(object, type = c("working", "standardised"), ...)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{type}{character string. The type of residual to be returned.}
  \item{...}{additional arguments.}
}

\value{
  Raw or standardised residuals extracted from 'object'.
}

\note{
  The 'standardised' residuals are available for least squares estimation
  with or without Box-Cox transformation or variance as a power of the
  mean.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())

## Displaying the residual plot (raw residuals)
plot(fitted(ryegrass.m1), residuals(ryegrass.m1))

## Using the standardised residuals
plot(fitted(ryegrass.m1), residuals(ryegrass.m1, type = "standard"))

}
\keyword{models}
\keyword{nonlinear}
