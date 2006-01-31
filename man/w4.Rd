\name{w4}

\alias{w4}

\title{The four-parameter Weibull function}

\description{
  'w4' provides the four-parameter Weibull function, self starter function and names of the parameters.  
}

\usage{
w4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), useDer = FALSE)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}

\details{
  The four-parameter Weibull model is given by the expression
  \deqn{ f(x) = c + (d-c)\exp(-\exp(b(\log(x)-e))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
}

\value{
  See \code{\link{weibull}}.
}

\references{Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{w3}}. A more flexible function, allowing fixing or constraining parameters, is \code{\link{weibull}}.}

\examples{

## Fitting a four-parameter Weibull model
model1 <- multdrc(rgr~dose, pct, data = TM, fct = w4())
summary(model1)


rm(model1)
}
\keyword{models}
\keyword{nonlinear}
