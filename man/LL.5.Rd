\name{LL.5}

\alias{LL.5}
\alias{l5}

\title{The five-parameter log-logistic function}

\description{
  'LL.5' provides the five-parameter log-logistic function, self starter function and names of the parameters.
}

\usage{
  LL.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)

  l5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{Additional arguments (see \code{\link{llogistic}}).}
}

\details{
  The five-parameter logistic function is given by the expression
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
  
  The function is asymmetric for \eqn{f} different from 1.
  
}
\value{
  See \code{\link{llogistic}}.
}
\references{
  Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, 
  \emph{Int. Statist. Rev.}, \bold{47}, 1--12.
}

\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Related functions are \code{\link{LL.4}} and \code{\link{LL.3}}.}
\examples{

model1 <- multdrc(ryegrass, fct = LL.5())
summary(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
