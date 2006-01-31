\name{ml4b}
\alias{ml4b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The modified logistic function}
\description{
  'ml4b' provides Cedergreen-Ritz-Streibig's modified logistic dose response function.
}
\usage{
ml4b(names = c("b", "c", "d", "e", "f"), useDer = FALSE)
}
\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable (see above).}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}
\details{
  The model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+f \exp(-1/(x^{0.5}))}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account.
}
\value{
  See \code{\link{mlogistic}}.
}
\references{Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved empirical models describing hormesis, Submitted.}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Similar functions are \code{\link{ml4a}} and \code{\link{ml4c}}.}
\examples{

## Modified logistic model
model1 <- multdrc(hormesis[,c(2,1)], fct=ml4b())
summary(model1)
ED(model1, c(50))

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
