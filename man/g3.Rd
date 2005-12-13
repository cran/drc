\name{g3}
\alias{g3}
\title{The three-parameter Gompertz function}
\description{
  'g3' provides the three-parameter Gompertz function, self starter function and names of the parameters.  
}
\usage{
g3(names = c("b", "d", "e"), useDer = FALSE)
}
\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}
\details{
  The three-parameter Gompertz model is given by the expression
  \deqn{ f(x) = 0 + (d-0)\exp(-\exp(b(\log(x)-e))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
}
\value{
  See \code{\link{gompertz}}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Related functions are \code{\link{g4}} and \code{\link{gompertz}}.}
\examples{

model1 <- multdrc(FA, fct=g3())
model1

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
