\name{logistic}

\alias{logistic}

\alias{L.3}
\alias{L.4}
\alias{L.5}

\alias{boltzmann}

\alias{B.3}
\alias{B.4}
\alias{B.5}
\alias{b3}
\alias{b4}
\alias{b5}

\title{The logistic (Boltzmann) model}

\description{
  The asymmetric logistic (Boltzmann) model for dose-response relationships.
}

\usage{
  logistic(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), fctName, fctText) 

  L.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
  L.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
  L.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))  

  boltzmann(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), fctName, fctText)

  B.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
  B.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
  B.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))  
}


\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}    
}

\details{
  The default arguments yields the five-parameter logistic mean function given by the expression
  
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}}
  
  The model is different from the log-logistic models \code{\link{llogistic}} and \code{\link{llogistic2}}
  where the term \deqn{log(x)} is used instead of \deqn{x}.
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

% \references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

## Fitting the logistic model
m1 <- multdrc(ryegrass, fct = logistic())
summary(m1)

## Fitting an asymmetric logistic model
##  requires installing the package 'NISTnls'
# m2 <- drm(y~x, data = Ratkowsky3, 
# fct = b5(fixed = c(NA, 0, NA, NA, NA)))
# plot(m2)
# summary(m2)  
## okay agreement with NIST values
##  for the two parameters that are the same

}
\keyword{models}
\keyword{nonlinear}
