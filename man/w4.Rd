\name{W1.4}

\alias{W1.4}
\alias{w4}
\alias{W2.4}

\title{The four-parameter Weibull functions}

\description{
  'W1.4' and 'W2.4' provide the four-parameter Weibull functions, self starter function and 
  names of the parameters.   
}

\usage{
  W1.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
  
  W2.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The equations for the mean functions are given at \code{\link{weibull1}}.
}

\value{
  See \code{\link{weibull1}}.
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).
}

\author{Christian Ritz}

\note{This function is for use with the model fitting function \code{\link{drm}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{W1.3}}. A more flexible function, allowing 
fixing or constraining parameters, is available through \code{\link{weibull1}}.}

\examples{

## Fitting a four-parameter Weibull (type 1) model
terbuthylazin.m1 <- drm(rgr~dose, data = terbuthylazin, fct = W1.4())
summary(terbuthylazin.m1)

}
\keyword{models}
\keyword{nonlinear}
