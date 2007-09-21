\name{rdrm}

\alias{rdrm}

\title{Simulating a dose-response curve}

\description{
  Simulation of a dose-response curve with user-specified dose values and error distribution.
}

\usage{
rdrm(nosim, fct, mpar, xerror, xpar = 1, yerror = "rnorm", ypar = c(0, 1), onlyY = FALSE, seedVal = 20070527)
}

\arguments{
  \item{nosim}{numeric. The number of simulated curves to be returned.}
  \item{fct}{list. Any built-in function in \code{\link{drc}} or a list with similar components.}
  \item{mpar}{numeric. The model parameters to be supplied to \code{fct}.}
  \item{xerror}{numeric or character. The distribution for the dose values.}
  \item{xpar}{numeric vector supplying the parameter values defining the distribution for the dose values. 
  If \code{xerror} is a distribution then remember that the number of dose values also is part of this argument 
  (the first argument).}
  \item{yerror}{numeric or character. The error distribution for the response values.}
  \item{ypar}{numeric vector supplying the parameter values defining the error distribution for the 
  response values.}
  \item{onlyY}{logical. If TRUE then only the response values are returned (useful in simulations). 
  Otherwise both dose values and response values (and for binomial data also the weights) are returned.}
  \item{seedVal}{numeric. The seed to be used.}
}

\details{
  The distribution for the dose values can either be a fixed set of dose values (a numeric vector)
  used repeatedly for creating all curves or be a distribution specified as a character string resulting in
  varying dose values from curve to curve.
  
  The error distribution for the response values can be any continuous distribution
  like \code{\link{rnorm}} or \code{\link{rgamma}}. Alternatively it can be the binomial distribution.
}

\value{
  A list with up to 3 components (depending on the value of the \code{onlyY} argument).
}

\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

## Model fit
m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())

## 10 random dose-response curves based on the model fit
sim10 <- rdrm(10, LL.4(), coef(m1), xerror = ryegrass$conc)
sim10

}

\keyword{models}
\keyword{nonlinear}
