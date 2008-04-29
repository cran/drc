\name{expDecay}

\alias{expDecay}

\title{Exponential decay model}

\description{
  Exponential decay model with a lower limit
}

\usage{
expDecay(fixed = c(NA, NA, NA), names = c("init", "plateau", "k"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The default parameter names are: init, plateau, k.}
}

\details{
  The exponential decay model is a three-parameter model with mean function:
  
  \deqn{ f(x) = init + (plateau-init)(\exp(-k x))}
  
  The parameter init is the upper limit (attained at x=0), the parameter plateau is the lower limit 
  reached for x going to infinity and the parameter k (>0) is determining the steepness of the 
  decay. Thus the curve is monotonously decreasing in x.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{}

\author{Christian Ritz}

%\note{}

\seealso{Similar models giving exponential increasing curves are \code{\link{AR.2}} and \code{\link{AR.3}}.}

\examples{

## Fitting an exponential decay model
ryegrass.m1<-drm(rootl~conc, data=ryegrass, fct=expDecay())

plot(ryegrass.m1)

summary(ryegrass.m1)

}
\keyword{models}
\keyword{nonlinear}
