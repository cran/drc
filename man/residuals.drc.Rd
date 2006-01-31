\name{residuals.drc}
\alias{residuals.drc}
\title{Extracting residuals from model}
\description{
  'residuals' extracts residuals from an object of class 'drc'.
}
\usage{
residuals.drc(object, ...)
}
\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}
\value{
  Rae residuals extracted from 'object'.
}
\author{Christian Ritz}
\examples{

model1 <- multdrc(FA)
plot(fitted(model1), residuals(model1))  # a residual plot

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
