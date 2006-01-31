\name{l2}

\alias{l2}

\title{The two-parameter logistic function}

\description{
  'l2' provides the two-parameter logistic function where the lower limit is fixed at 0 and the upper limit
  is fixed at 1, mostly suitable for binomial/quantal responses. 
}

\usage{
  l2(fixed = c(NA, NA), names = c("b", "e"), useDer = FALSE)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}

\details{
  The two-parameter logistic function is given by the expression
  \deqn{ f(x) = \frac{1}{1+\exp(b(\log(x)-\log(e)))}}
  
  The function is symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{logistic}}.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{multdrc}}.
}

\seealso{
  Related functions are \code{\link{l3}}, \code{\link{l4}}, \code{\link{l5}} and the more general 
  \code{\link{logistic}}.
}

\examples{

## Fitting a two-parameter logistic model 
##  to binomial responses (a logit model)
model1 <- multdrc(number/total~dose, weights=total, 
data=earthworms, fct = l2(), type = "binomial")
plot(model1)  # not fitting at the upper limit!

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
