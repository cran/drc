\name{l4}
\alias{l4}

\title{The four-parameter logistic function}

\description{
  'l4' provides the four-parameter logistic function, self starter function, names of the parameters and, optionally,
  first and second derivatives for a faster estimation.  
}

\usage{
l4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), useDer = FALSE, w = FALSE)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, if FALSE (default) they are not supplied.}
  \item{w}{logical. If TRUE weighted regression is applied to determine the start values.}  
}

\details{
  The four-parameter logistic function is given by the expression
  \deqn{ f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}}
  
  The function is symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{logistic}}.
}

\references{
   Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (p. 330).
}

\author{Christian Ritz and Jens C. Streibig}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{l3}}. See also \code{\link{l5}}.}

\examples{

model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)
model1

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
