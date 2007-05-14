\name{llogistic}

\alias{llogistic}

\title{The log-logistic function}

\description{
  'llogistic' provides a very general way of specifying log-logistic models,
  under various constraints on the parameters.
}

\usage{
  llogistic(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf),
  fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"),
  ss = c("1", "2", "3"), ssfct = NULL)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
    NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f 
    (see under 'Details').}
  \item{ss}{character string indicating the self starter function to use.}
  \item{ssfct}{a self starter function to be used.}  
}

\details{
  The default arguments yields the five-parameter log-logistic function given by the expression
  
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
  
  If the parameter \eqn{f} differs from 1 then the function is asymmetric; otherwise it
  is symmetric (on log scale).
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

\references{
  Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, 
  \emph{Int. Statist. Rev.}, \bold{47}, 1--12.

  Seber, G. A. F. and Wild, C. J. (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (p. 330).
}

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{drm}}.
}

\seealso{
  For convenience several special cases are available: 
  \code{\link{LL.3}}, \code{\link{LL.4}} and \code{\link{LL.5}}. 
}

\examples{


}
\keyword{models}
\keyword{nonlinear}
