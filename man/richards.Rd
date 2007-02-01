\name{richards}

\alias{richards}

\title{Richards model}

\description{
  'richards' 
}

\usage{
  richards(fixed = c(NA, NA, NA, NA, NA), names = c("c", "d", "delta", "kappa", "gamma"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: c, d, delta, kappa, gamma (see under 'Details').}}
}

\details{
  The model is given by the expression
  
  \deqn{ f(x) = c + (d-c) ( 1 + (delta-1)\exp(-kappa(\log(x)-\log(gamma))) )^(1/(1-delta))}
  
  It is a five-parameter model suitable for asymmetric dose-response relationships. In fact
  it is a reparameterisation of the five-parameter log-logistic model.
  
  For delta=2 the model simplifies the four-parameter log-logistic model.
}

\value{
  An object containing with all components needed in \code{\link{multdrc}}.
}

\references{
  Seber, G. A. F. and Wild, C. J. (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 332--337).
}

\author{Christian Ritz}

\note{
  It may be difficult to get this model to converge.
}

%\seealso{}

\examples{

## Fitting Richards model
m1 <- multdrc(ryegrass, fct = richards())

}
\keyword{models}
\keyword{nonlinear}
