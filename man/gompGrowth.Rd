\name{gompGrowth}

\alias{gompGrowth.1}
\alias{gompGrowth.2}
\alias{gompGrowth.3}

\title{Gompertz growth models}

\description{
  Gompertz growth model, with biologically meaningful parameters. Different parameterisations have been
  included for specific cases and needs.
}

\usage{
gompGrowth.1(fixed = c(NA, NA, NA), names = c("c", "m", "plateau"))
gompGrowth.2(fixed = c(NA, NA, NA), names = c("c", "d", "plateau"))
gompGrowth.3(fixed = c(NA, NA, NA), names = c("b", "c", "plateau"))

}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The default parameter names are: init, m, plateau.}
}

\details{
  The Gompertz growth model is a Gompertz curve, that has been reparameterised to include some biologically 
  meaningful parameters. The mean function for gompGrowth.1() is:
  
  \deqn{ f(x) = f(x) = plateau * exp ( - (m/c) * exp ( - c *  x ) )}
  
  The parameter plateau is the final plant weight, reached for x going to infinity
  the parameter c is relative growth rate at inflection point and
  the parameter m is the initial relative growth rate (when x=0).
  Thus the curve is monotonously increasing in x.
  The mean function for gompGrowth.2() is:

  \deqn{ f(x) = plateau * exp ( - exp ( c * ( d - x ) )) }
  
  where the parameter c is the relative growth rate at inflection point and the parameter d is the abscissa of the 
  inflection point.
  The mean function for gompGrowth.3() is the classical Gompertz function:
  
  \deqn{ f(x) = plateau * exp ( - b * exp ( - c *  x ) ) }
  
  where b is proportional to the initial relative growth rate (m = b * c).
  
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names.
}

\references{Roderick Hunt, 1982. Plant Growth Curves. Edward Arnold Publisher, Great Britain, 248 pp}

\author{Andrea Onofri}

\note{Growth functions are generally fitted on log-transformed weight data, which equals to setting bc parameter to 0}

%\seealso{}

\examples{

## Fitting a Gompertz growth curve

beet.model <- drm(weightInf ~ DAE, data  = beetGrowth, fct=gompGrowth.1())
plot(beet.model, log = "")
summary(beet.model)

}
\keyword{models}
\keyword{nonlinear}

\concept{growth Gompertz}
