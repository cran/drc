\name{fitted.drc}
\alias{fitted.drc}
\title{Extract fitted values from model}
\description{
  Extracts fitted values from an object of class 'drc'. 
}
\usage{
fitted.drc(object, ...)
}
\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}
\value{
  Fitted values extracted from 'object'.
}
\author{Christian Ritz}

\examples{

model1 <- multdrc(FA)
plot(fitted(model1), residuals(model1))  # a residual plot

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
