\name{boltzmann}

\alias{boltzmann}

\alias{B.3}
\alias{B.4}
\alias{B.5}
\alias{b3}
\alias{b4}
\alias{b5}

\title{The Boltzmann model}

\description{
  The Boltzmann model for dose-response relationships.
}

\usage{
  boltzmann(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))

  B.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))

  B.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
  
  B.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))  
}


\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b, c, d, e, f (see under 'Details').}}
}

\details{
  The model is also called the logistic model (different from the log-logistic model!).
}

\value{
  An object containing with all components needed in \code{\link{multdrc}}.
}

% \references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

## Fitting the Boltzmann model
m1 <- multdrc(ryegrass, fct = boltzmann())
summary(m1)

## Fitting an asymmetric Boltzmann model
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
