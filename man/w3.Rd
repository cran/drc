\name{w3}

\alias{w3}

\title{The three-parameter Weibull function}

\description{
  'w3' provides the three-parameter Weibull function, self starter function and names of the parameters.  
}

\usage{
w3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), useDer = FALSE)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}

\details{
  The three-parameter Weibull model is given by the expression
  \deqn{ f(x) = 0 + (d-0)\exp(-\exp(b(\log(x)-e))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
}

\value{
  See \code{\link{weibull}}.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Related functions are \code{\link{w4}} and \code{\link{weibull}}.}

\examples{

## Fitting a three-parameter Weibull model
model1 <- multdrc(FA, fct = w3())
model1

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
