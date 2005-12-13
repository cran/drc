\name{summary.drc}
\alias{summary.drc}
\title{Summarising non-linear model fits}
\description{
  'summary' compiles a comprehensive summary for objects of class 'drc'. 
}
\usage{
summary.drc(object, od = FALSE, ...)
}
\arguments{
  \item{object}{an object of class 'drc'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}  
  \item{...}{additional arguments.}
}

\value{
  A list of summary statistics, including parameter estimates and estimated standard errors, of the fitted 'drc' model given by \code{object}.
}
\author{Christian Ritz}

\examples{

model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)
summary(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
