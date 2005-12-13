\name{l3}
\alias{l3}
\title{The three-parameter logistic function (lower limit equal to zero)}
\description{
  'l3' provides the three-parameter logistic function, self starter function and names of the parameters. 
}
\usage{
l3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), useDer = FALSE)
}
\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}
\details{
  The three-parameter logistic function is given by the expression
  \deqn{ f(x) = 0 + \frac{d-0}{1+\exp(b(\log(x)-\log(e)))}}
  
  The function is symmetric about the inflection point (\eqn{e}).  

}
\value{
  See \code{\link{logistic}}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Christian}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{A related function is \code{\link{l4}} and the more general \code{\link{logistic}}.}
\examples{

model1 <- multdrc(FA, fct=l3())
model1

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
