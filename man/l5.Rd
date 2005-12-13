\name{l5}
\alias{l5}
\title{The five-parameter logistic function}
\description{
  'l5' provides the five-parameter logistic function, self starter function and names of the parameters.
}
\usage{
l5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), useDer = FALSE)
}
\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}
\details{
  The five-parameter logistic function is given by the expression
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
  
  The function is asymmetric for \eqn{f} different from 1.
  
}
\value{
  See \code{\link{logistic}}.
}
\references{Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, \emph{Int. Statist. Rev.}, \bold{47}, 1--12.}

\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Related functions are \code{\link{l4}} and \code{\link{l3}}.}
\examples{

model1 <- multdrc(FA, fct=l5())
summary(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
